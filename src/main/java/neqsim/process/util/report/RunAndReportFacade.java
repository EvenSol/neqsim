package neqsim.process.util.report;

import java.nio.charset.StandardCharsets;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Base64;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.UUID;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import neqsim.process.equipment.ProcessEquipmentInterface;
import neqsim.process.processmodel.ProcessSystem;
import neqsim.process.util.report.RunAndReportRequest.RunMode;
import neqsim.thermo.system.SystemInterface;

/**
 * Facade that runs a flowsheet and emits a predictable KPI + provenance payload.
 */
public class RunAndReportFacade {
  private static final Logger logger = LogManager.getLogger(RunAndReportFacade.class);

  /**
   * Execute a steady-state or transient run, capture KPIs and provenance, and serialize a report.
   *
   * @param system flowsheet to execute
   * @param request configuration for the run
   * @return standardized result payload
   */
  public RunReportResult runAndReport(ProcessSystem system, RunAndReportRequest request) {
    if (system == null) {
      throw new IllegalArgumentException("Process system cannot be null");
    }
    RunAndReportRequest effectiveRequest = request == null ? new RunAndReportRequest() : request;

    String flowsheetHash = computeFlowsheetHash(system, effectiveRequest.getHashConfig());
    UUID runId = effectiveRequest.getRunId();

    if (effectiveRequest.getRunMode() == RunMode.TRANSIENT) {
      executeTransient(system, effectiveRequest, runId);
    } else {
      system.run(runId);
    }

    Map<String, Object> kpis = collectKpis(system);
    RunReportResult.Provenance provenance = new RunReportResult.Provenance(flowsheetHash,
        effectiveRequest.getTemplateId(), collectThermoModelVersions(system));

    String serializedReport = new Report(system).generateJsonReport(effectiveRequest.getReportConfig());

    return new RunReportResult(effectiveRequest.getRunMode(), runId, kpis, provenance,
        serializedReport);
  }

  private void executeTransient(ProcessSystem system, RunAndReportRequest request, UUID runId) {
    for (int step = 0; step < request.getTransientSteps(); step++) {
      system.runTransient(request.getTimeStep(), runId);
    }
  }

  private Map<String, Object> collectKpis(ProcessSystem system) {
    Map<String, Object> kpis = new LinkedHashMap<>();

    kpis.put("unitOperationCount", system.getUnitOperations().size());
    kpis.put("power_kW", system.getPower("kW"));
    kpis.put("heaterDuty_kW", system.getHeaterDuty("kW"));
    kpis.put("coolerDuty_kW", system.getCoolerDuty("kW"));

    double exergyChange = system.getExergyChange("J");
    kpis.put("exergyChange_kJ", exergyChange / 1.0e3);
    kpis.put("maxMassBalanceError_percent", calculateMaxMassBalanceError(system));

    return kpis;
  }

  private double calculateMaxMassBalanceError(ProcessSystem system) {
    double maxError = Double.NaN;
    try {
      for (ProcessSystem.MassBalanceResult result : system.checkMassBalance().values()) {
        double percentError = result.getPercentError();
        if (Double.isNaN(percentError)) {
          continue;
        }
        if (Double.isNaN(maxError) || Math.abs(percentError) > maxError) {
          maxError = Math.abs(percentError);
        }
      }
    } catch (Exception ex) {
      logger.warn("Unable to calculate mass balance KPI", ex);
    }
    return maxError;
  }

  private List<String> collectThermoModelVersions(ProcessSystem system) {
    Set<String> versions = new LinkedHashSet<>();
    for (ProcessEquipmentInterface unit : system.getUnitOperations()) {
      try {
        SystemInterface thermo = unit.getThermoSystem();
        if (thermo != null) {
          String modelName = thermo.getModelName();
          String implementationVersion = thermo.getClass().getPackage().getImplementationVersion();
          if (implementationVersion != null && !implementationVersion.isEmpty()) {
            versions.add(modelName + "@" + implementationVersion);
          } else {
            versions.add(modelName);
          }
        }
      } catch (Exception ex) {
        logger.debug("Unable to read thermo model from unit {}", unit.getName(), ex);
      }
    }
    return new ArrayList<>(versions);
  }

  private String computeFlowsheetHash(ProcessSystem system, ReportConfig hashConfig) {
    String report = new Report(system).generateJsonReport(hashConfig);
    try {
      MessageDigest digest = MessageDigest.getInstance("SHA-256");
      byte[] hash = digest.digest(report.getBytes(StandardCharsets.UTF_8));
      return Base64.getEncoder().encodeToString(hash);
    } catch (NoSuchAlgorithmException e) {
      throw new IllegalStateException("SHA-256 algorithm not available", e);
    }
  }
}
