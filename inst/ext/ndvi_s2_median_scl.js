//VERSION=3
// Sentinel-2 L2A NDVI: per-pixel median over all clear observations in the
// requested time range.
// Clear = scene classification (SCL) 4 vegetation, 5 not vegetated, 6 water.
// Excluded: 0 no data, 1 saturated, 2 dark/topographic shadow, 3 cloud shadow,
// 7 unclassified, 8/9 cloud medium/high probability, 10 thin cirrus, 11 snow.
// Output: band 1 median NDVI (NaN if no clear observation),
//         band 2 number of clear observations.
function setup() {
  return {
    input: [{
      bands: ["B04", "B08", "SCL", "dataMask"],
      units: ["REFLECTANCE", "REFLECTANCE", "DN", "DN"]
    }],
    output: { id: "default", bands: 2, sampleType: "FLOAT32" },
    mosaicking: "ORBIT"
  };
}

const CLEAR = [4, 5, 6];

function median(v) {
  v.sort(function (a, b) { return a - b; });
  const m = Math.floor(v.length / 2);
  return v.length % 2 ? v[m] : (v[m - 1] + v[m]) / 2;
}

function evaluatePixel(samples) {
  const v = [];
  for (let i = 0; i < samples.length; i++) {
    const s = samples[i];
    if (s.dataMask !== 1 || CLEAR.indexOf(s.SCL) < 0) continue;
    const d = s.B08 + s.B04;
    if (!(d > 0)) continue;
    v.push((s.B08 - s.B04) / d);
  }
  if (!v.length) return [NaN, 0];
  return [median(v), v.length];
}
