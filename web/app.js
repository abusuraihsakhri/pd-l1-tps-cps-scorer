/* Local-only PD-L1 arithmetic. No network requests or persistent storage. */
(function (root) {
  "use strict";
  function count(value, name, positive) {
    var n = Number(value);
    if (String(value).trim() === "" || !Number.isSafeInteger(n) || n < (positive ? 1 : 0)) {
      throw new Error(name + " must be a " + (positive ? "positive" : "non-negative") + " whole number.");
    }
    return n;
  }
  function area(value, name, positive) {
    var n = Number(value);
    if (String(value).trim() === "" || !Number.isFinite(n) || n < 0 || (positive && n === 0)) {
      throw new Error(name + " must be a " + (positive ? "positive" : "non-negative") + " finite number.");
    }
    return n;
  }
  function compute(values) {
    var pos = count(values.positive, "PD-L1-positive tumor cells", false);
    var total = count(values.total, "Total viable tumor cells", true);
    var lymph = count(values.lymphocytes, "Positive lymphocytes", false);
    var macro = count(values.macrophages, "Positive macrophages", false);
    if (pos > total) throw new Error("Positive tumor cells cannot exceed total viable tumor cells.");
    var tpsRaw = 100 * pos / total;
    var cpsRaw = 100 * (pos + lymph + macro) / total;
    var cps = Math.min(100, cpsRaw);
    return {
      tps: tpsRaw, tpsCategory: tpsRaw < 1 ? "Below 1%" : tpsRaw < 50 ? "1% to <50%" : "50% or higher",
      cps: cps, cpsCategory: cps < 1 ? "Below 1" : cps < 10 ? "1 to <10" : cps < 50 ? "10 to <50" : "50 or higher",
      capped: cpsRaw > 100,
      lowCount: total < 100,
      totalPositive: pos + lymph + macro
    };
  }
  function computeIC(positiveArea, tumorArea) {
    var pos = area(positiveArea, "PD-L1-positive immune-cell area", false);
    var total = area(tumorArea, "Total tumor area", true);
    if (pos > total) throw new Error("Positive immune-cell area cannot exceed total tumor area.");
    var score = 100 * pos / total;
    return { ic: score, category: score < 1 ? "IC0" : score < 5 ? "IC1" : "IC2/3" };
  }
  // CSV parser handles quoted commas, escaped quotes and CRLF.
  function parseCSV(input) {
    var rows = [], row = [], field = "", quoted = false;
    input = input.replace(/^\uFEFF/, "");
    for (var i = 0; i < input.length; i++) {
      var ch = input[i];
      if (quoted) {
        if (ch === '"' && input[i + 1] === '"') { field += '"'; i++; }
        else if (ch === '"') quoted = false;
        else field += ch;
      } else if (ch === '"') {
        if (field !== "") throw new Error("Invalid CSV quoting.");
        quoted = true;
      } else if (ch === ",") { row.push(field); field = ""; }
      else if (ch === "\n" || ch === "\r") {
        if (ch === "\r" && input[i + 1] === "\n") i++;
        row.push(field); field = "";
        if (row.some(function (s) { return s.trim() !== ""; })) rows.push(row);
        row = [];
      } else { field += ch; }
    }
    if (quoted) throw new Error("Unterminated quoted CSV field.");
    row.push(field);
    if (row.some(function (s) { return s.trim() !== ""; })) rows.push(row);
    return rows;
  }
  function csvEscape(s) {
    var val = String(s);
    // Prevent spreadsheet formula interpretation of any exported text fields.
    if (/^[\s]*[=+@-]/.test(val)) val = "'" + val;
    return '"' + val.replace(/"/g, '""') + '"';
  }
  function batchFromCSV(content, limit) {
    var rows = parseCSV(content);
    if (!rows.length) throw new Error("The CSV is empty.");
    var required = ["pd_l1_positive_tumor_cells", "total_viable_tumor_cells"];
    var headers = rows[0].map(function (h) { return h.trim(); });
    required.forEach(function (h) { if (!headers.includes(h)) throw new Error("Missing required column: " + h); });
    var maxRows = limit || 1000;
    if (rows.length - 1 > maxRows) throw new Error("File exceeds the " + maxRows + "-case limit.");
    var results = [];
    rows.slice(1).forEach(function (cols, i) {
      if (cols.length !== headers.length) throw new Error("CSV row " + (i + 2) + " has an unexpected column count.");
      var obj = {};
      headers.forEach(function (h, j) { obj[h] = cols[j].trim(); });
      try {
        var v = compute({
          positive: obj.pd_l1_positive_tumor_cells,
          total: obj.total_viable_tumor_cells,
          lymphocytes: obj.pd_l1_positive_lymphocytes === undefined || obj.pd_l1_positive_lymphocytes === "" ? "0" : obj.pd_l1_positive_lymphocytes,
          macrophages: obj.pd_l1_positive_macrophages === undefined || obj.pd_l1_positive_macrophages === "" ? "0" : obj.pd_l1_positive_macrophages
        });
        results.push({ row: i + 2, tps: v.tps, cps: v.cps, capped: v.capped, adequacy: v.lowCount ? "Check assay minimum viable tumor-cell count" : "" });
      } catch (error) { throw new Error("CSV row " + (i + 2) + ": " + error.message); }
    });
    return results;
  }
  function exportCSV(rows) {
    var lines = [["source_row", "tps_percent", "cps", "cps_capped", "adequacy_note"].map(csvEscape).join(",")];
    rows.forEach(function (r) {
      lines.push([r.row, r.tps.toFixed(2), r.cps.toFixed(2), r.capped ? "yes" : "no", r.adequacy].map(csvEscape).join(","));
    });
    return lines.join("\r\n") + "\r\n";
  }
  var scorer = { compute: compute, computeIC: computeIC, parseCSV: parseCSV, batchFromCSV: batchFromCSV, exportCSV: exportCSV };
  if (typeof module !== "undefined" && module.exports) module.exports = scorer;
  root.PDL1Scorer = scorer;

  if (typeof document === "undefined") return;
  document.addEventListener("DOMContentLoaded", function () {
    var form = document.getElementById("score-form");
    var out = document.getElementById("score-result");
    var icForm = document.getElementById("ic-form");
    var icOut = document.getElementById("ic-result");
    var batchFile = document.getElementById("batch-file");
    var batchStatus = document.getElementById("batch-status");
    var download = document.getElementById("download-results");
    var batchRows = [];
    var byId = function (id) { return document.getElementById(id); };
    var show = function (el, message, fail) {
      el.textContent = message; el.classList.toggle("error", !!fail);
    };
    var fmt = function (n) { return n.toFixed(2); };
    form.addEventListener("submit", function (event) {
      event.preventDefault();
      try {
        var r = compute({
          positive: byId("tumor-pos").value, total: byId("tumor-total").value,
          lymphocytes: byId("lymph-pos").value, macrophages: byId("macro-pos").value
        });
        show(out, "TPS: " + fmt(r.tps) + "% (" + r.tpsCategory + ")\n" +
          "CPS: " + fmt(r.cps) + " (" + r.cpsCategory + ")\n" +
          "PD-L1-positive cells in CPS numerator: " + r.totalPositive +
          (r.capped ? "\nCPS capped at 100 per the 22C3 scoring definition." : "") +
          (r.lowCount ? "\nAdequacy warning: fewer than 100 viable tumor cells. Check assay-specific adequacy criteria." : ""));
      } catch (error) { show(out, error.message, true); }
    });
    icForm.addEventListener("submit", function (event) {
      event.preventDefault();
      try {
        var r = computeIC(byId("ic-positive").value, byId("ic-total").value);
        show(icOut, "IC: " + fmt(r.ic) + "% (" + r.category + ")\nArea-based estimate; apply assay-specific interpretation.");
      } catch (error) { show(icOut, error.message, true); }
    });
    batchFile.addEventListener("change", async function () {
      batchRows = []; download.disabled = true;
      var file = batchFile.files && batchFile.files[0];
      if (!file) { show(batchStatus, "Select a CSV file."); return; }
      if (file.size > 2 * 1024 * 1024) { show(batchStatus, "Maximum file size is 2 MB.", true); return; }
      try {
        batchRows = batchFromCSV(await file.text(), 1000);
        show(batchStatus, "Validated " + batchRows.length + " case(s). Results are processed locally. " +
          (batchRows.length ? "Use Download results to save TPS and CPS." : "No data rows."));
        download.disabled = batchRows.length === 0;
      } catch (error) { show(batchStatus, error.message, true); }
    });
    download.addEventListener("click", function () {
      if (!batchRows.length) return;
      var blob = new Blob([exportCSV(batchRows)], { type: "text/csv;charset=utf-8" });
      var url = URL.createObjectURL(blob);
      var a = document.createElement("a");
      a.href = url; a.download = "pd-l1-scores.csv"; document.body.appendChild(a); a.click(); a.remove();
      setTimeout(function () { URL.revokeObjectURL(url); }, 0);
    });
  });
})(typeof window !== "undefined" ? window : globalThis);
