//VERSION=3
function setup() {
  return {
    input: [{ bands: ["B04","B08","dataMask"], units: "REFLECTANCE" }],
    output: { id: "default", bands: 2, sampleType: "FLOAT32" } // default!
  };
}
// function evaluatePixel(s) {
//  const ndvi = (s.B08 - s.B04) / (s.B08 + s.B04);
//  return [ndvi, s.dataMask];  // Band1: NDVI (float), Band2: Mask (0/1 als float)
// }

function evaluatePixel(s) {
  const ndvi = (s.B08 - s.B04) / (s.B08 + s.B04);
  const dm = s.dataMask === 1;
  return [dm ? ndvi : NaN, s.dataMask]; // NDVI band carries NaN where invalid
}

