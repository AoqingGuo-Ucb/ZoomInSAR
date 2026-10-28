const $ = (id) => document.getElementById(id);
let projects = [], selected, currentMap, map, velocityOverlay, projectOutline, overviewLayer;
const chartCanvas = $("chart");
const californiaBounds = [[32.45, -124.55], [42.15, -114.05]];

function color(value, limit) {
  if (value === null || Number.isNaN(value)) return [0, 0, 0, 0];
  const intensity = Math.min(1, Math.abs(value) / limit);
  const end = value < 0 ? [35, 71, 162] : [191, 38, 52];
  return end.map((channel) => Math.round(245 + (channel - 245) * intensity)).concat(225);
}
function showStatus(text) { $("status").textContent = text; }
async function api(path) { const response = await fetch(path); const body = await response.json(); if (!response.ok) throw new Error(body.error || "Request failed"); return body; }

function initializeMap() {
  map = L.map("map", { zoomControl: true, minZoom: 5, maxZoom: 18, maxBounds: [[30, -127], [44, -111]] });
  L.tileLayer("https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}", {
    attribution: "Tiles &copy; Esri — Source: Esri, Maxar, Earthstar Geographics, and the GIS User Community", maxZoom: 18,
  }).addTo(map);
  map.fitBounds(californiaBounds);
  map.on("click", clickMap);
}

function overlayImage(payload, enlargeSparsePixels = false) {
  const canvas = document.createElement("canvas"); canvas.width = payload.width; canvas.height = payload.height;
  const context = canvas.getContext("2d"); const image = context.createImageData(payload.width, payload.height);
  payload.values.forEach((value, index) => image.data.set(color(value, payload.limit_mm_year), index * 4));
  context.putImageData(image, 0, 0);
  const validCount = payload.values.reduce((count, value) => count + (value === null ? 0 : 1), 0);
  if (enlargeSparsePixels && validCount > 0 && validCount < 120000) {
    const radius = validCount < 30000 ? 3 : 2;
    payload.values.forEach((value, index) => {
      if (value === null) return;
      const [red, green, blue] = color(value, payload.limit_mm_year);
      context.fillStyle = `rgb(${red}, ${green}, ${blue})`;
      context.beginPath();
      context.arc(index % payload.width + 0.5, Math.floor(index / payload.width) + 0.5, radius, 0, 2 * Math.PI);
      context.fill();
      context.strokeStyle = "rgba(25, 35, 45, 0.45)";
      context.lineWidth = 0.6;
      context.stroke();
    });
  }
  return canvas.toDataURL("image/png");
}

function drawMap(payload) {
  currentMap = payload;
  const [west, south, east, north] = payload.bounds, bounds = [[south, west], [north, east]];
  if (overviewLayer) { map.removeLayer(overviewLayer); overviewLayer = undefined; }
  if (velocityOverlay) map.removeLayer(velocityOverlay);
  if (projectOutline) map.removeLayer(projectOutline);
  velocityOverlay = L.imageOverlay(overlayImage(payload, true), bounds, { opacity: 0.84, interactive: false, className: "velocity-overlay" }).addTo(map);
  projectOutline = L.rectangle(bounds, { color: "#f4cd39", weight: 2, fill: false, interactive: false }).addTo(map);
  $("limit-left").textContent = payload.limit_mm_year.toFixed(1);
  $("limit-right").textContent = payload.limit_mm_year.toFixed(1);
  $("period").textContent = `${payload.period[0]} — ${payload.period[1]}`;
}

function drawOverview(payload) {
  currentMap = undefined;
  if (velocityOverlay) { map.removeLayer(velocityOverlay); velocityOverlay = undefined; }
  if (projectOutline) { map.removeLayer(projectOutline); projectOutline = undefined; }
  if (overviewLayer) map.removeLayer(overviewLayer);
  overviewLayer = L.layerGroup();
  payload.projects.forEach((project) => {
    const [west, south, east, north] = project.bounds, bounds = [[south, west], [north, east]];
    L.imageOverlay(overlayImage(project), bounds, { opacity: 0.72, interactive: false }).addTo(overviewLayer);
    L.rectangle(bounds, { color: "#f4cd39", weight: 1, fill: false, interactive: false }).bindTooltip(project.name).addTo(overviewLayer);
    L.circleMarker([(south + north) / 2, (west + east) / 2], {
      radius: 11, color: "#ffffff", weight: 3, fillColor: "#e22f3b", fillOpacity: 0.98,
    }).bindTooltip(project.name, { direction: "top", offset: [0, -10] }).on("click", () => {
      $("project").value = project.id;
      chooseProject();
    }).addTo(overviewLayer);
  });
  overviewLayer.addTo(map);
  $("map-title").textContent = "All California results";
  $("period").textContent = "Full available time spans";
  $("limit-left").textContent = payload.limit_mm_year.toFixed(1);
  $("limit-right").textContent = payload.limit_mm_year.toFixed(1);
  map.fitBounds(californiaBounds);
}

function zoomToProject() {
  if (!currentMap) return;
  const [west, south, east, north] = currentMap.bounds;
  map.fitBounds([[south, west], [north, east]], { padding: [30, 30], maxZoom: 14 });
}

async function loadMap() {
  if (!selected) return;
  showStatus("Calculating and mapping velocity…");
  try {
    drawMap(await api(`/api/velocity?project=${encodeURIComponent(selected.id)}&start=${$("start").value}&end=${$("end").value}&resolution=720`));
    showStatus("Velocity overlay loaded. Zoom using the map controls; click inside the project outline for a time series.");
  } catch (error) { showStatus(`Error: ${error.message}`); }
}

async function loadOverview() {
  showStatus("Loading all California velocity results…");
  try {
    drawOverview(await api("/api/overview"));
    showStatus("All available results are shown. Select a project to inspect its time range and point time series.");
  } catch (error) { showStatus(`Error: ${error.message}`); }
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
  if (!currentMap) { showStatus("Select a project before requesting a point time series."); return; }
  const [west, south, east, north] = currentMap.bounds, { lat, lng } = event.latlng;
  if (lng < west || lng > east || lat < south || lat > north) { showStatus("Click inside the selected project outline to retrieve a point time series."); return; }
  showStatus("Loading point time series…");
  try { const series = await api(`/api/timeseries?project=${encodeURIComponent(selected.id)}&lon=${lng}&lat=${lat}`); drawChart(series); $("series-title").textContent = "Cumulative LOS displacement"; $("point").textContent = `${series.longitude.toFixed(5)}, ${series.latitude.toFixed(5)}`; showStatus(`Nearest radar pixel: row ${series.row}, column ${series.col}`); }
  catch (error) { showStatus(`Error: ${error.message}`); }
}

function configureProject() {
  selected = projects.find((project) => project.id === $("project").value); if (!selected) return;
  $("start").value = selected.dates[0]; $("end").value = selected.dates.at(-1);
  $("details").innerHTML = `<h2>Project details</h2><p>${selected.shape[1]} × ${selected.shape[0]} pixels<br>${selected.dates.length} acquisitions<br>${selected.dates[0]} to ${selected.dates.at(-1)}</p>`;
}

function chooseProject() {
  configureProject();
  if (!selected) return;
  $("map-title").textContent = selected.name;
  loadMap();
}

async function initialize() {
  initializeMap();
  try {
    projects = (await api("/api/projects")).projects;
    const select = $("project"); projects.forEach((project) => select.add(new Option(project.name, project.id)));
    select.onchange = chooseProject; $("update").onclick = loadMap; $("zoom-project").onclick = zoomToProject;
    $("show-overview").onclick = loadOverview;
    $("reset").onclick = () => { $("start").value = selected.dates[0]; $("end").value = selected.dates.at(-1); loadMap(); };
    configureProject();
    loadOverview();
  } catch (error) { showStatus(`Could not load dashboard: ${error.message}`); }
}
initialize();
