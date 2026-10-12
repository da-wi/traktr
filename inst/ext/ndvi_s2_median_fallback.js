//VERSION=3
// Sentinel-2 L2A NDVI with a tiered composite over the requested time range.
// Per pixel, the first tier that has observations is used:
//   1 clear:   median NDVI of SCL 4 vegetation, 5 not vegetated, 6 water
//   2 relaxed: median NDVI of tier-1 classes plus SCL 2 dark area,
//              7 unclassified, 10 thin cirrus (still excludes dense cloud,
//              cloud shadow and snow)
//   3 maximum: maximum NDVI of all observations with data (maximum value
//              composite; clouds and haze lower NDVI, so the maximum is the
//              least contaminated observation)
// Output: band 1 NDVI (NaN only if there is no observation at all),
//         band 2 number of observations in the tier used,
//         band 3 tier (1, 2, 3; NaN if no observation).
function setup() {
  return {
    input: [{
      bands: ["B04", "B08", "SCL", "dataMask"],
      units: ["REFLECTANCE", "REFLECTANCE", "DN", "DN"]
    }],
    output: { id: "default", bands: 3, sampleType: "FLOAT32" },
    mosaicking: "ORBIT"
  };
}

const CLEAR   = [4, 5, 6];
const RELAXED = [2, 4, 5, 6, 7, 10];

function median(v) {
  v.sort(function (a, b) { return a - b; });
  const m = Math.floor(v.length / 2);
  return v.length % 2 ? v[m] : (v[m - 1] + v[m]) / 2;
}

function evaluatePixel(samples) {
  const clear = [], relaxed = [];
  let max = -Infinity, n = 0;
  for (let i = 0; i < samples.length; i++) {
    const s = samples[i];
    if (s.dataMask !== 1) continue;
    const d = s.B08 + s.B04;
    if (!(d > 0)) continue;
    const v = (s.B08 - s.B04) / d;
    n++;
    if (v > max) max = v;
    if (RELAXED.indexOf(s.SCL) >= 0) {
      relaxed.push(v);
      if (CLEAR.indexOf(s.SCL) >= 0) clear.push(v);
    }
  }
  if (clear.length)   return [median(clear), clear.length, 1];
  if (relaxed.length) return [median(relaxed), relaxed.length, 2];
  if (n)              return [max, n, 3];
  return [NaN, 0, NaN];
}
