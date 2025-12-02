package neqsim.process.util.report;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonObject;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.UUID;

/**
 * Standardized payload from the run-and-report facade.
 */
public class RunReportResult {
  /**
   * Lightweight provenance record that can be persisted together with KPI results.
   */
  public static class Provenance {
    private final String flowsheetHash;
    private final String templateId;
    private final List<String> thermoModelVersions;

    public Provenance(String flowsheetHash, String templateId, List<String> thermoModelVersions) {
      this.flowsheetHash = flowsheetHash;
      this.templateId = templateId;
      this.thermoModelVersions = thermoModelVersions;
    }

    public String getFlowsheetHash() {
      return flowsheetHash;
    }

    public String getTemplateId() {
      return templateId;
    }

    public List<String> getThermoModelVersions() {
      return thermoModelVersions;
    }
  }

  private final RunAndReportRequest.RunMode runMode;
  private final UUID runId;
  private final Map<String, Object> kpis;
  private final Provenance provenance;
  private final String report;

  public RunReportResult(RunAndReportRequest.RunMode runMode, UUID runId, Map<String, Object> kpis,
      Provenance provenance, String report) {
    this.runMode = runMode;
    this.runId = runId;
    this.kpis = new LinkedHashMap<>(kpis);
    this.provenance = provenance;
    this.report = report;
  }

  public RunAndReportRequest.RunMode getRunMode() {
    return runMode;
  }

  public UUID getRunId() {
    return runId;
  }

  public Map<String, Object> getKpis() {
    return Collections.unmodifiableMap(kpis);
  }

  public Provenance getProvenance() {
    return provenance;
  }

  public String getReport() {
    return report;
  }

  /**
   * Convert the result into a JSON string using a predictable schema.
   *
   * @return JSON representation of the run result
   */
  public String toJson() {
    Gson gson = new GsonBuilder().serializeSpecialFloatingPointValues().setPrettyPrinting().create();
    return gson.toJson(toJsonObject());
  }

  /**
   * Convert the result to a JSON object for further processing.
   *
   * @return structured JSON object
   */
  public JsonObject toJsonObject() {
    JsonObject obj = new JsonObject();
    obj.addProperty("runMode", runMode.name());
    obj.addProperty("runId", runId.toString());

    JsonObject kpiObject = new JsonObject();
    for (Map.Entry<String, Object> entry : kpis.entrySet()) {
      Object value = entry.getValue();
      if (value instanceof Number) {
        kpiObject.addProperty(entry.getKey(), (Number) value);
      } else if (value instanceof Boolean) {
        kpiObject.addProperty(entry.getKey(), (Boolean) value);
      } else {
        kpiObject.addProperty(entry.getKey(), value == null ? null : value.toString());
      }
    }
    obj.add("kpis", kpiObject);

    if (provenance != null) {
      JsonObject prov = new JsonObject();
      prov.addProperty("flowsheetHash", provenance.getFlowsheetHash());
      prov.addProperty("templateId", provenance.getTemplateId());

      JsonObject thermo = new JsonObject();
      int counter = 0;
      for (String version : provenance.getThermoModelVersions()) {
        thermo.addProperty("model" + counter++, version);
      }
      obj.add("provenance", prov);
      prov.add("thermoModels", thermo);
    }

    obj.addProperty("report", report);
    return obj;
  }
}
