"use strict";
const test = require("node:test");
const assert = require("node:assert/strict");
const scorer = require("./app.js");
test("TPS and CPS calculations and cap", () => {
  assert.equal(scorer.compute({ positive: 50, total: 100, lymphocytes: 5, macrophages: 3 }).tps, 50);
  assert.equal(scorer.compute({ positive: 50, total: 100, lymphocytes: 5, macrophages: 3 }).cps, 58);
  const capped = scorer.compute({ positive: 100, total: 100, lymphocytes: 30, macrophages: 5 });
  assert.equal(capped.cps, 100); assert.equal(capped.capped, true);
});
test("thresholds use unrounded score", () => {
  const result = scorer.compute({ positive: 12, total: 1250, lymphocytes: 0, macrophages: 0 });
  assert.equal(result.tpsCategory, "Below 1%");
  assert.equal(result.cpsCategory, "Below 1");
  assert.equal(scorer.compute({ positive: 50, total: 100, lymphocytes: 0, macrophages: 0 }).tpsCategory, "50% or higher");
});
test("reject invalid count and area combinations", () => {
  assert.throws(() => scorer.compute({ positive: 101, total: 100, lymphocytes: 0, macrophages: 0 }));
  assert.throws(() => scorer.compute({ positive: "1.5", total: 100, lymphocytes: 0, macrophages: 0 }));
  assert.throws(() => scorer.compute({ positive: 1, total: 0, lymphocytes: 0, macrophages: 0 }));
  assert.throws(() => scorer.computeIC("Infinity", 100));
  assert.throws(() => scorer.computeIC(101, 100));
  assert.deepEqual(scorer.computeIC(5, 100), { ic: 5, category: "IC2/3" });
});
test("CSV quoted cells, validation and export", () => {
  assert.deepEqual(scorer.parseCSV('a,b\r\n"quoted, value","escaped ""quote"""\r\n'), [["a", "b"], ["quoted, value", 'escaped "quote"']]);
  const rows = scorer.batchFromCSV("pd_l1_positive_tumor_cells,total_viable_tumor_cells\r\n0,100\r\n", 1000);
  assert.equal(rows.length, 1); assert.equal(rows[0].cps, 0);
  assert.match(scorer.exportCSV(rows), /tps_percent/);
  assert.throws(() => scorer.batchFromCSV("bad,column\n1,2\n"), /Missing required column/);
  assert.throws(() => scorer.batchFromCSV("pd_l1_positive_tumor_cells,total_viable_tumor_cells\n1,0\n"), /CSV row 2/);
});
