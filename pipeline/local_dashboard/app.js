const $ = (id) => document.getElementById(id);
let projects = [], profiles = {}, selected, currentMap, map, opticalMap, satelliteLayer, velocityOverlay, projectOutline, opticalOutline, overviewLayer, pointMarker, opticalPointMarker, synchronizing = false, selectedPoint, plottedSeries;
const chartCanvas = $("chart");
const californiaBounds = [[32.45, -124.55], [42.15, -114.05]];

function color(value, limit, mode = "diverging") {
  if (value === null || Number.isNaN(value)) return [0, 0, 0, 0];
  if (mode === "quality") {
    const strength = Math.min(1, Math.max(0, value / limit));
    return [Math.round(55 + 170 * strength), Math.round(45 + 165 * strength), Math.round(100 - 45 * strength), 225];
  }
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
  const limit = payload.limit_mm_year ?? payload.limit;
  payload.values.forEach((value, index) => image.data.set(color(value, limit, payload.color_mode), index * 4));
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
    const [red, green, blue] = color(value, limit, payload.color_mode);
    context.fillStyle = `rgb(${red}, ${green}, ${blue})`;
    context.strokeStyle = "rgba(25, 35, 45, 0.42)";
    context.lineWidth = 0.45;
    context.beginPath();
    context.arc(index % payload.width + 0.5, Math.floor(index / payload.width) + 0.5, nearestRadius(index), 0, 2 * Math.PI);
    context.fill(); context.stroke();
  });
  return canvas.toDataURL("image/png");
}

function clearSelectedLayers(keepSelectedPoint = false) {
  if (velocityOverlay) { map.removeLayer(velocityOverlay); velocityOverlay = undefined; }
  if (projectOutline) { map.removeLayer(projectOutline); projectOutline = undefined; }
  if (opticalOutline && opticalMap) { opticalMap.removeLayer(opticalOutline); opticalOutline = undefined; }
  if (!keepSelectedPoint) {
    if (pointMarker) { map.removeLayer(pointMarker); pointMarker = undefined; }
    if (opticalPointMarker && opticalMap) { opticalMap.removeLayer(opticalPointMarker); opticalPointMarker = undefined; }
  }
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
  clearSelectedLayers(Boolean(selectedPoint));
  const bounds = boundsFor(payload);
  velocityOverlay = L.imageOverlay(makeVelocityImage(payload, true), bounds, { opacity: 0.90, interactive: false, className: "velocity-overlay" }).addTo(map);
  projectOutline = L.rectangle(bounds, { color: "#f4cd39", weight: 2, fill: false, interactive: false }).addTo(map);
  opticalOutline = L.rectangle(bounds, { color: "#f4cd39", weight: 2, fill: false, interactive: false }).addTo(opticalMap);
  const pot = $("product").value === "pot";
  const limit = payload.limit_mm_year ?? payload.limit;
  $("map-title").textContent = pot ? `${selected.name} — Pixel-offset tracking` : selected.name;
  $("limit-left").textContent = pot && payload.color_mode === "quality" ? "0" : `-${limit.toFixed(1)} ${payload.unit || "mm/year"}`;
  $("limit-right").textContent = pot && payload.color_mode === "quality" ? "1.0 correlation" : `+${limit.toFixed(1)} ${payload.unit || "mm/year"}`;
  $("period").textContent = `${payload.period[0]} — ${payload.period[1]}`;
  if (!pot && selectedPoint) showSelectedPoint(selectedPoint.latitude, selectedPoint.longitude);
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
  $("limit-left").textContent = `-${payload.limit_mm_year.toFixed(1)} mm/year`;
  $("limit-right").textContent = `+${payload.limit_mm_year.toFixed(1)} mm/year`;
  requestAnimationFrame(() => map.fitBounds(californiaBounds));
}

function zoomToProject() { if (currentMap) fitProject(currentMap); }
async function loadMap(zoomToProject = false) {
  if (!selected) return;
  showStatus("Calculating and mapping velocity…");
  try {
    drawMap(await api(`/api/velocity?project=${encodeURIComponent(selected.id)}&start=${$("start").value}&end=${$("end").value}&resolution=720`), zoomToProject);
    if (selectedPoint) {
      await loadPointTimeSeries(selectedPoint.longitude, selectedPoint.latitude);
      showStatus("Velocity and the selected point time series were updated for the chosen period.");
    } else {
      showStatus("Velocity and optical maps are synchronized. Click either map inside the yellow boundary for a time series.");
    }
  }
  catch (error) { showStatus(`Error: ${error.message}`); }
}
async function loadPotMap(zoomToProject = false) {
  if (!selected?.pixel_offset) return;
  showStatus("Loading pixel-offset result…");
  selectedPoint = undefined;
  try {
    const pair = $("pot-pair").value || "0", field = $("pot-field").value;
    const payload = await api(`/api/pixel-offset?project=${encodeURIComponent(selected.id)}&pair=${pair}&field=${field}&resolution=720`);
    drawMap(payload, zoomToProject);
    $("series-title").textContent = "Pixel-offset tracking result";
    $("point").textContent = "POT is pairwise; it does not provide a cumulative LOS time series.";
    chartCanvas.getContext("2d").clearRect(0, 0, chartCanvas.width, chartCanvas.height);
    $("chart-readout").textContent = "Switch to Phase InSAR to inspect a point displacement time series.";
    $("series-values").querySelector("div").innerHTML = "";
    showStatus(`POT ${field} map loaded for ${payload.period[0]} to ${payload.period[1]}.`);
  } catch (error) { showStatus(`Error: ${error.message}`); }
}
async function loadOverview() {
  showStatus("Loading all California velocity results…");
  try { drawOverview(await api("/api/overview")); showStatus("All available results are shown. Click a red project marker to open its velocity and optical maps."); }
  catch (error) { showStatus(`Error: ${error.message}`); }
}

function renderSeriesValues(series) {
  const rows = series.dates.map((date, index) => {
    const value = series.displacement_mm[index];
    return `<tr><td>${date}</td><td>${value === null ? "No data" : `${value.toFixed(1)} mm`}</td></tr>`;
  }).join("");
  $("series-values").querySelector("div").innerHTML = `<table><thead><tr><th>Date</th><th>Cumulative LOS displacement</th></tr></thead><tbody>${rows}</tbody></table>`;
}

function drawChart(series) {
  const valid = series.displacement_mm.map((value, index) => value === null ? null : [index, value]).filter(Boolean);
  plottedSeries = { series, valid, pad: 45, width: 900, height: 260 };
  chartCanvas.width = 900; chartCanvas.height = 260;
  const context = chartCanvas.getContext("2d"), width = chartCanvas.width, height = chartCanvas.height, pad = 45;
  context.clearRect(0, 0, width, height); if (!valid.length) return;
  const values = valid.map((point) => point[1]), ymin = Math.min(...values), ymax = Math.max(...values), span = Math.max(ymax - ymin, 1);
  context.strokeStyle = "#d8dddb"; context.lineWidth = 1;
  for (let i = 0; i < 5; i += 1) { const y = pad + i * (height - 2 * pad) / 4; context.beginPath(); context.moveTo(pad, y); context.lineTo(width - pad, y); context.stroke(); context.fillStyle = "#60717e"; context.font = "12px system-ui"; context.fillText((ymax - i * span / 4).toFixed(1), 4, y + 4); }
  context.strokeStyle = "#b83c3c"; context.lineWidth = 2; context.beginPath();
  valid.forEach(([index, value], number) => { const x = pad + index * (width - 2 * pad) / Math.max(series.dates.length - 1, 1), y = height - pad - (value - ymin) * (height - 2 * pad) / span; number ? context.lineTo(x, y) : context.moveTo(x, y); }); context.stroke();
  context.fillStyle = "#b83c3c";
  valid.forEach(([index, value]) => { const x = pad + index * (width - 2 * pad) / Math.max(series.dates.length - 1, 1), y = height - pad - (value - ymin) * (height - 2 * pad) / span; context.beginPath(); context.arc(x, y, 2.5, 0, 2 * Math.PI); context.fill(); });
  context.fillStyle = "#60717e"; context.font = "12px system-ui"; context.fillText(series.dates[0], pad, height - 12); context.fillText(series.dates.at(-1), width - pad - 78, height - 12); context.fillText("Cumulative LOS displacement (mm)", pad, 16);
  renderSeriesValues(series);
}

function inspectChart(event) {
  if (!plottedSeries?.valid.length) return;
  const rect = chartCanvas.getBoundingClientRect();
  const x = (event.clientX - rect.left) * chartCanvas.width / rect.width;
  const { series, valid, pad, width } = plottedSeries;
  const index = Math.round((x - pad) * Math.max(series.dates.length - 1, 1) / (width - 2 * pad));
  const closest = valid.reduce((best, point) => Math.abs(point[0] - index) < Math.abs(best[0] - index) ? point : best);
  if (Math.abs(closest[0] - index) > 2) return;
  $("chart-readout").textContent = `${series.dates[closest[0]]}: ${closest[1].toFixed(1)} mm cumulative LOS displacement`;
}

async function loadPointTimeSeries(longitude, latitude) {
  const series = await api(`/api/timeseries?project=${encodeURIComponent(selected.id)}&lon=${longitude}&lat=${latitude}&start=${$("start").value}&end=${$("end").value}`);
  selectedPoint = { latitude: series.latitude, longitude: series.longitude };
  showSelectedPoint(series.latitude, series.longitude);
  drawChart(series);
  $("series-title").textContent = `Cumulative LOS displacement: ${series.dates[0]} — ${series.dates.at(-1)}`;
  $("point").textContent = `${series.longitude.toFixed(5)}, ${series.latitude.toFixed(5)}`;
  return series;
}

async function clickMap(event) {
  if ($("product").value === "pot") { showStatus("POT is a pairwise offset product. Switch to Phase InSAR to request a point time series."); return; }
  if (!currentMap) { showStatus("Click a red project marker or choose a project before requesting a point time series."); return; }
  const [[south, west], [north, east]] = boundsFor(currentMap), { lat, lng } = event.latlng;
  if (lng < west || lng > east || lat < south || lat > north) { showStatus("Click inside the yellow project boundary to retrieve a point time series."); return; }
  showStatus("Loading point time series…");
  try { const series = await loadPointTimeSeries(lng, lat); showStatus(`Nearest radar pixel: row ${series.row}, column ${series.col}`); }
  catch (error) { showStatus(`Error: ${error.message}`); }
}

function configureProject(id) {
  if (selected?.id !== id) {
    selectedPoint = undefined;
    plottedSeries = undefined;
    chartCanvas.getContext("2d").clearRect(0, 0, chartCanvas.width, chartCanvas.height);
    $("series-title").textContent = "Click the velocity map";
    $("point").textContent = "";
    $("chart-readout").textContent = "Hover over a plotted observation to inspect its date and displacement.";
    $("series-values").querySelector("div").innerHTML = "";
  }
  selected = projects.find((project) => project.id === id); if (!selected) return false;
  $("project").value = selected.id;
  const product = $("product"); product.innerHTML = "";
  product.add(new Option("Phase InSAR", "insar"));
  if (selected.pixel_offset) product.add(new Option("Pixel-offset tracking (POT)", "pot"));
  product.value = "insar";
  const pot = selected.pixel_offset;
  $("pot-pair").innerHTML = "";
  if (pot) pot.pairs.forEach((pair) => $("pot-pair").add(new Option(`${pair.label} (${pair.baseline_days} days)`, pair.index)));
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
function loadSelectedProduct(zoomToProject = false) {
  $("phase-controls").hidden = $("product").value === "pot";
  $("pot-controls").hidden = $("product").value !== "pot";
  if ($("product").value === "pot") loadPotMap(zoomToProject); else loadMap(zoomToProject);
}
function selectProject(id, zoomToProject) { if (configureProject(id)) loadSelectedProduct(zoomToProject); }

async function initialize() {
  initializeMap();
  try {
    const [projectResponse, profileResponse] = await Promise.all([api("/api/projects"), api("/api/profiles")]);
    projects = projectResponse.projects; profiles = profileResponse.profiles;
    const select = $("project"); projects.forEach((project) => select.add(new Option(project.name, project.id)));
    select.onchange = () => selectProject(select.value, true);
    $("update").onclick = () => loadMap(false); $("update-pot").onclick = () => loadPotMap(false); $("zoom-project").onclick = zoomToProject; $("show-overview").onclick = loadOverview;
    $("product").onchange = () => loadSelectedProduct(true);
    $("reset").onclick = () => { if (!selected) return; $("start").value = selected.dates[0]; $("end").value = selected.dates.at(-1); loadMap(false); };
    chartCanvas.addEventListener("mousemove", inspectChart);
    configureProject(projects[0].id); loadOverview();
  } catch (error) { showStatus(`Could not load dashboard: ${error.message}`); }
}
initialize();
