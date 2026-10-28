const $ = (id) => document.getElementById(id);
let projects = [], profiles = {}, selected, currentMap, map, opticalMap, satelliteLayer, velocityOverlay, projectOutline, opticalOutline, overviewLayer, pointMarker, opticalPointMarker, synchronizing = false;
const chartCanvas = $("chart");
const californiaBounds = [[32.45, -124.55], [42.15, -114.05]];

function color(value, limit) {
  if (value === null || Number.isNaN(value)) return [0, 0, 0, 0];
  const strength = Math.min(1, Math.abs(value) / limit);
  const endpoint = value < 0 ? [35, 71, 162] : [191, 38, 52];
  return endpoint.map((channel) => Math.round(245 + (channel - 245) * strength)).concat(225);
}
function showStatus(text) { $("status").textContent = text; }
async function api(path) { const response = await fetch(path); const body = await response.json(); if (!response.ok) throw new Error(body.error || "Request failed"); return body; }
function imageryLayer() { return L.tileLayer("https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}", { attribution: "Tiles &copy; Esri — Source: Esri, Maxar, Earthstar Geographics, and the GIS User Community", maxZoom: 18 }); }
function boundsFor(payload) { const [west, south, east, north] = payload.bounds; return [[south, west], [north, east]]; }

function initializeMap() {
  map = L.map("map", { zoomControl: true, minZoom: 5, maxZoom: 18, maxBounds: [[30, -127], [44, -111]] });
  satelliteLayer = imageryLayer().addTo(map);
  map.fitBounds(californiaBounds);
  map.on("click", clickMap);
}

function initializeOpticalMap() {
  if (opticalMap) return;
  opticalMap = L.map("optical-map", { zoomControl: true, minZoom: 5, maxZoom: 18, maxBounds: [[30, -127], [44, -111]] });
  imageryLayer().addTo(opticalMap);
  opticalMap.on("click", clickMap);
  const synchronize = (source, target) => () => {
    if (synchronizing || $("map-layout").classList.contains("overview")) return;
    synchronizing = true;
    target.setView(source.getCenter(), source.getZoom(), { animate: false });
    synchronizing = false;
  };
  map.on("moveend", synchronize(map, opticalMap));
  opticalMap.on("moveend", synchronize(opticalMap, map));
}

function enterProjectLayout() {
  $("map-layout").classList.remove("overview");
  initializeOpticalMap();
  requestAnimationFrame(() => { map.invalidateSize(); opticalMap.invalidateSize(); });
}

function enterOverviewLayout() {
  $("map-layout").classList.add("overview");
  if (!map.hasLayer(satelliteLayer)) satelliteLayer.addTo(map);
  requestAnimationFrame(() => map.invalidateSize());
}

function makeVelocityImage(payload, drawDots = false) {
  const canvas = document.createElement("canvas"); canvas.width = payload.width; canvas.height = payload.height;
  const context = canvas.getContext("2d"); const image = context.createImageData(payload.width, payload.height);
  payload.values.forEach((value, index) => image.data.set(color(value, payload.limit_mm_year), index * 4));
  context.putImageData(image, 0, 0);
  if (!drawDots) return canvas.toDataURL("image/png");

  // Enlarge isolated valid pixels only. The nearest occupied cell limits each
  // radius, so adjacent velocity observations never overlap.
  const occupied = new Uint8Array(payload.values.length);
  payload.values.forEach((value, index) => { if (value !== null) occupied[index] = 1; });
  const nearestRadius = (index) => {
    const row = Math.floor(index / payload.width), column = index % payload.width;
    let nearestSquared = Infinity;
    for (let dy = -6; dy <= 6; dy += 1) for (let dx = -6; dx <= 6; dx += 1) {
      if (dx === 0 && dy === 0) continue;
      const y = row + dy, x = column + dx;
      if (y < 0 || y >= payload.height || x < 0 || x >= payload.width || !occupied[y * payload.width + x]) continue;
      nearestSquared = Math.min(nearestSquared, dx * dx + dy * dy);
    }
    return Math.min(3, Math.max(0.43, Math.sqrt(nearestSquared) / 2 - 0.08));
  };
  payload.values.forEach((value, index) => {
    if (value === null) return;
    const [red, green, blue] = color(value, payload.limit_mm_year);
    context.fillStyle = `rgb(${red}, ${green}, ${blue})`;
    context.strokeStyle = "rgba(25, 35, 45, 0.42)";
    context.lineWidth = 0.45;
    context.beginPath();
    context.arc(index % payload.width + 0.5, Math.floor(index / payload.width) + 0.5, nearestRadius(index), 0, 2 * Math.PI);
    context.fill(); context.stroke();
  });
  return canvas.toDataURL("image/png");
}

function clearSelectedLayers() {
  if (velocityOverlay) { map.removeLayer(velocityOverlay); velocityOverlay = undefined; }
  if (projectOutline) { map.removeLayer(projectOutline); projectOutline = undefined; }
  if (opticalOutline && opticalMap) { opticalMap.removeLayer(opticalOutline); opticalOutline = undefined; }
  if (pointMarker) { map.removeLayer(pointMarker); pointMarker = undefined; }
  if (opticalPointMarker && opticalMap) { opticalMap.removeLayer(opticalPointMarker); opticalPointMarker = undefined; }
}

function showSelectedPoint(latitude, longitude) {
  if (pointMarker) map.removeLayer(pointMarker);
  if (opticalPointMarker && opticalMap) opticalMap.removeLayer(opticalPointMarker);
  const icon = () => L.divIcon({ className: "series-arrow", html: "▼", iconSize: [30, 34], iconAnchor: [15, 28] });
  pointMarker = L.marker([latitude, longitude], { icon: icon(), interactive: false }).addTo(map);
  opticalPointMarker = L.marker([latitude, longitude], { icon: icon(), interactive: false }).addTo(opticalMap);
}

function fitProject(payload) {
  const bounds = boundsFor(payload);
  map.fitBounds(bounds, { padding: [28, 28], maxZoom: 15 });
  opticalMap.fitBounds(bounds, { padding: [28, 28], maxZoom: 15 });
}

function drawMap(payload, zoomToProject) {
  currentMap = payload;
  enterProjectLayout();
  if (overviewLayer) { map.removeLayer(overviewLayer); overviewLayer = undefined; }
  clearSelectedLayers();
  const bounds = boundsFor(payload);
  velocityOverlay = L.imageOverlay(makeVelocityImage(payload, true), bounds, { opacity: 0.90, interactive: false, className: "velocity-overlay" }).addTo(map);
  projectOutline = L.rectangle(bounds, { color: "#f4cd39", weight: 2, fill: false, interactive: false }).addTo(map);
  opticalOutline = L.rectangle(bounds, { color: "#f4cd39", weight: 2, fill: false, interactive: false }).addTo(opticalMap);
  $("map-title").textContent = selected.name;
  $("limit-left").textContent = payload.limit_mm_year.toFixed(1);
  $("limit-right").textContent = payload.limit_mm_year.toFixed(1);
  $("period").textContent = `${payload.period[0]} — ${payload.period[1]}`;
  if (zoomToProject) requestAnimationFrame(() => fitProject(payload));
}

function drawOverview(payload) {
  currentMap = undefined;
  enterOverviewLayout(); clearSelectedLayers();
  if (overviewLayer) map.removeLayer(overviewLayer);
  overviewLayer = L.layerGroup();
  payload.projects.forEach((project) => {
    const bounds = boundsFor(project);
    L.imageOverlay(makeVelocityImage(project), bounds, { opacity: 0.72, interactive: false }).addTo(overviewLayer);
    L.rectangle(bounds, { color: "#f4cd39", weight: 1, fill: false, interactive: false }).bindTooltip(project.name).addTo(overviewLayer);
    L.circleMarker([(bounds[0][0] + bounds[1][0]) / 2, (bounds[0][1] + bounds[1][1]) / 2], { radius: 11, color: "#ffffff", weight: 3, fillColor: "#e22f3b", fillOpacity: 0.98 }).bindTooltip(project.name, { direction: "top", offset: [0, -10] }).on("click", () => selectProject(project.id, true)).addTo(overviewLayer);
  });
  overviewLayer.addTo(map);
  $("map-title").textContent = "All California results";
  $("period").textContent = "Full available time spans";
  $("limit-left").textContent = payload.limit_mm_year.toFixed(1);
  $("limit-right").textContent = payload.limit_mm_year.toFixed(1);
  requestAnimationFrame(() => map.fitBounds(californiaBounds));
}

function zoomToProject() { if (currentMap) fitProject(currentMap); }
async function loadMap(zoomToProject = false) {
  if (!selected) return;
  showStatus("Calculating and mapping velocity…");
  try { drawMap(await api(`/api/velocity?project=${encodeURIComponent(selected.id)}&start=${$("start").value}&end=${$("end").value}&resolution=720`), zoomToProject); showStatus("Velocity and optical maps are synchronized. Click either map inside the yellow boundary for a time series."); }
  catch (error) { showStatus(`Error: ${error.message}`); }
}
async function loadOverview() {
  showStatus("Loading all California velocity results…");
  try { drawOverview(await api("/api/overview")); showStatus("All available results are shown. Click a red project marker to open its velocity and optical maps."); }
  catch (error) { showStatus(`Error: ${error.message}`); }
}

function drawChart(series) {
  const valid = series.displacement_mm.map((value, index) => value === null ? null : [index, value]).filter(Boolean);
  chartCanvas.width = 900; chartCanvas.height = 260;
  const context = chartCanvas.getContext("2d"), width = chartCanvas.width, height = chartCanvas.height, pad = 45;
  context.clearRect(0, 0, width, height); if (!valid.length) return;
  const values = valid.map((point) => point[1]), ymin = Math.min(...values), ymax = Math.max(...values), span = Math.max(ymax - ymin, 1);
  context.strokeStyle = "#d8dddb"; context.lineWidth = 1;
  for (let i = 0; i < 5; i += 1) { const y = pad + i * (height - 2 * pad) / 4; context.beginPath(); context.moveTo(pad, y); context.lineTo(width - pad, y); context.stroke(); context.fillStyle = "#60717e"; context.font = "12px system-ui"; context.fillText((ymax - i * span / 4).toFixed(1), 4, y + 4); }
  context.strokeStyle = "#b83c3c"; context.lineWidth = 2; context.beginPath();
  valid.forEach(([index, value], number) => { const x = pad + index * (width - 2 * pad) / (series.dates.length - 1), y = height - pad - (value - ymin) * (height - 2 * pad) / span; number ? context.lineTo(x, y) : context.moveTo(x, y); }); context.stroke();
  context.fillStyle = "#60717e"; context.font = "12px system-ui"; context.fillText(series.dates[0], pad, height - 12); context.fillText(series.dates.at(-1), width - pad - 78, height - 12); context.fillText("Cumulative LOS displacement (mm)", pad, 16);
}

async function clickMap(event) {
  if (!currentMap) { showStatus("Click a red project marker or choose a project before requesting a point time series."); return; }
  const [[south, west], [north, east]] = boundsFor(currentMap), { lat, lng } = event.latlng;
  if (lng < west || lng > east || lat < south || lat > north) { showStatus("Click inside the yellow project boundary to retrieve a point time series."); return; }
  showStatus("Loading point time series…");
  try { const series = await api(`/api/timeseries?project=${encodeURIComponent(selected.id)}&lon=${lng}&lat=${lat}`); showSelectedPoint(series.latitude, series.longitude); drawChart(series); $("series-title").textContent = "Cumulative LOS displacement"; $("point").textContent = `${series.longitude.toFixed(5)}, ${series.latitude.toFixed(5)}`; showStatus(`Nearest radar pixel: row ${series.row}, column ${series.col}`); }
  catch (error) { showStatus(`Error: ${error.message}`); }
}

function configureProject(id) {
  selected = projects.find((project) => project.id === id); if (!selected) return false;
  $("project").value = selected.id;
  $("start").value = selected.dates[0]; $("end").value = selected.dates.at(-1);
  $("details").innerHTML = `<h2>Project details</h2><p>${selected.shape[1]} × ${selected.shape[0]} pixels<br>${selected.dates.length} acquisitions<br>${selected.dates[0]} to ${selected.dates.at(-1)}</p>`;
  const profile = profiles[selected.id];
  if (!profile) {
    $("site-profile").innerHTML = "<h2>Landfill information</h2><p>No reviewed site profile has been loaded.</p>";
  } else {
    const facility = profile.facility, match = profile.match;
    const review = match.confidence === "high" ? "Coordinate match" : "Coordinate match — review before citation";
    const source = profile.sources?.[0];
    $("site-profile").innerHTML = `<h2>Landfill information</h2><p><strong>${facility.name || "Unresolved facility"}</strong><br>${facility.county || ""} County<br>Status: ${facility.operational_status || "unknown"}<br>Operator: ${facility.operator || "not listed"}<br>LEA: ${facility.enforcement_agency || "not listed"}<br><small>${review}; ${match.distance_km} km from project centre.</small>${source ? `<br><a href="${source.url}" target="_blank" rel="noreferrer">SWIS source</a>` : ""}</p>`;
  }
  return true;
}
function selectProject(id, zoomToProject) { if (configureProject(id)) loadMap(zoomToProject); }

async function initialize() {
  initializeMap();
  try {
    const [projectResponse, profileResponse] = await Promise.all([api("/api/projects"), api("/api/profiles")]);
    projects = projectResponse.projects; profiles = profileResponse.profiles;
    const select = $("project"); projects.forEach((project) => select.add(new Option(project.name, project.id)));
    select.onchange = () => selectProject(select.value, true);
    $("update").onclick = () => loadMap(false); $("zoom-project").onclick = zoomToProject; $("show-overview").onclick = loadOverview;
    $("reset").onclick = () => { if (!selected) return; $("start").value = selected.dates[0]; $("end").value = selected.dates.at(-1); loadMap(false); };
    configureProject(projects[0].id); loadOverview();
  } catch (error) { showStatus(`Could not load dashboard: ${error.message}`); }
}
initialize();
