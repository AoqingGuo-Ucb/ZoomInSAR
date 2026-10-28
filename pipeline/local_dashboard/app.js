const $ = (id) => document.getElementById(id);
let projects = [], selected, currentMap;
const mapCanvas = $("map"), chartCanvas = $("chart");

function color(value, limit) {
  if (value === null || Number.isNaN(value)) return [235, 232, 224, 0];
  const t = Math.min(1, Math.abs(value) / limit);
  const low = value < 0 ? [35, 71, 162] : [191, 38, 52];
  return low.map((channel) => Math.round(245 + (channel - 245) * t)).concat(225);
}
function showStatus(text) { $("status").textContent = text; }
async function api(path) { const r = await fetch(path); const body = await r.json(); if (!r.ok) throw new Error(body.error || "Request failed"); return body; }

function drawMap(payload) {
  currentMap = payload; const {width, height, values, limit_mm_year: limit} = payload;
  mapCanvas.width = width; mapCanvas.height = height; const ctx = mapCanvas.getContext("2d"), image = ctx.createImageData(width, height);
  values.forEach((v, i) => image.data.set(color(v, limit), i * 4)); ctx.putImageData(image, 0, 0);
  $("limit-left").textContent = limit.toFixed(1); $("limit-right").textContent = limit.toFixed(1);
  $("period").textContent = `${payload.period[0]} — ${payload.period[1]}`;
}
async function loadMap() {
  if (!selected) return; showStatus("Calculating and mapping velocity…");
  const start = $("start").value, end = $("end").value;
  try { drawMap(await api(`/api/velocity?project=${encodeURIComponent(selected.id)}&start=${start}&end=${end}`)); showStatus("Click the map to inspect a point time series."); }
  catch (e) { showStatus(`Error: ${e.message}`); }
}
function drawChart(series) {
  const valid = series.displacement_mm.map((v,i)=>v===null?null:[i,v]).filter(Boolean);
  chartCanvas.width = 900; chartCanvas.height = 260; const ctx = chartCanvas.getContext("2d"), w=chartCanvas.width,h=chartCanvas.height,p=45;
  ctx.clearRect(0,0,w,h); if (!valid.length) return;
  const ys=valid.map(x=>x[1]), ymin=Math.min(...ys), ymax=Math.max(...ys), span=Math.max(ymax-ymin,1);
  ctx.strokeStyle="#d8dddb"; ctx.lineWidth=1; for(let i=0;i<5;i++){const y=p+i*(h-2*p)/4;ctx.beginPath();ctx.moveTo(p,y);ctx.lineTo(w-p,y);ctx.stroke();ctx.fillStyle="#60717e";ctx.font="12px system-ui";ctx.fillText((ymax-i*span/4).toFixed(1),4,y+4)}
  ctx.strokeStyle="#b83c3c";ctx.lineWidth=2;ctx.beginPath(); valid.forEach(([i,v],n)=>{const x=p+i*(w-2*p)/(series.dates.length-1),y=h-p-(v-ymin)*(h-2*p)/span;n?ctx.lineTo(x,y):ctx.moveTo(x,y)});ctx.stroke();
  ctx.fillStyle="#60717e";ctx.font="12px system-ui";ctx.fillText(series.dates[0],p,h-12);ctx.fillText(series.dates.at(-1),w-p-78,h-12);ctx.fillText("Cumulative LOS displacement (mm)",p,16);
}
async function clickMap(event) {
  if (!currentMap) return; const rect=mapCanvas.getBoundingClientRect(), x=(event.clientX-rect.left)/rect.width, y=(event.clientY-rect.top)/rect.height;
  const [west,south,east,north]=currentMap.bounds, lon=west+x*(east-west), lat=north-y*(north-south);
  showStatus("Loading point time series…");
  try { const s=await api(`/api/timeseries?project=${encodeURIComponent(selected.id)}&lon=${lon}&lat=${lat}`); drawChart(s); $("series-title").textContent="Cumulative LOS displacement"; $("point").textContent=`${s.longitude.toFixed(5)}, ${s.latitude.toFixed(5)}`; showStatus(`Nearest radar pixel: row ${s.row}, column ${s.col}`); }
  catch(e){showStatus(`Error: ${e.message}`)}
}
function chooseProject() { selected=projects.find(p=>p.id===$("project").value); if(!selected)return; $("start").value=selected.dates[0];$("end").value=selected.dates.at(-1);$("map-title").textContent=selected.name;$("details").innerHTML=`<h2>Project details</h2><p>${selected.shape[1]} × ${selected.shape[0]} pixels<br>${selected.dates.length} acquisitions<br>${selected.dates[0]} to ${selected.dates.at(-1)}</p>`;loadMap(); }
async function initialize(){try{projects=(await api("/api/projects")).projects;const select=$("project");projects.forEach(p=>select.add(new Option(p.name,p.id)));select.onchange=chooseProject;$("update").onclick=loadMap;$("reset").onclick=()=>{ $("start").value=selected.dates[0];$("end").value=selected.dates.at(-1);loadMap();};mapCanvas.onclick=clickMap;chooseProject()}catch(e){showStatus(`Could not load dashboard: ${e.message}`)}}
initialize();
