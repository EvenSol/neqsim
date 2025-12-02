package neqsim.process.util.report;

import java.util.UUID;

/**
 * Configuration for the run-and-report facade.
 */
public class RunAndReportRequest {
  /** Supported run modes. */
  public enum RunMode {
    /** Execute a steady-state calculation. */
    STEADY_STATE,
    /** Execute a transient calculation for a configured number of steps. */
    TRANSIENT
  }

  private RunMode runMode = RunMode.STEADY_STATE;
  private int transientSteps = 1;
  private double timeStep = 1.0;
  private String templateId = "";
  private ReportConfig reportConfig = new ReportConfig(ReportConfig.DetailLevel.SUMMARY);
  private ReportConfig hashConfig = new ReportConfig(ReportConfig.DetailLevel.MINIMUM);
  private UUID runId = UUID.randomUUID();

  public RunMode getRunMode() {
    return runMode;
  }

  public RunAndReportRequest setRunMode(RunMode runMode) {
    this.runMode = runMode;
    return this;
  }

  public int getTransientSteps() {
    return transientSteps;
  }

  public RunAndReportRequest setTransientSteps(int transientSteps) {
    if (transientSteps < 1) {
      throw new IllegalArgumentException("Transient steps must be at least 1");
    }
    this.transientSteps = transientSteps;
    return this;
  }

  public double getTimeStep() {
    return timeStep;
  }

  public RunAndReportRequest setTimeStep(double timeStep) {
    if (timeStep <= 0.0) {
      throw new IllegalArgumentException("Time step must be positive");
    }
    this.timeStep = timeStep;
    return this;
  }

  public String getTemplateId() {
    return templateId;
  }

  public RunAndReportRequest setTemplateId(String templateId) {
    this.templateId = templateId == null ? "" : templateId;
    return this;
  }

  public ReportConfig getReportConfig() {
    return reportConfig;
  }

  public RunAndReportRequest setReportConfig(ReportConfig reportConfig) {
    if (reportConfig != null) {
      this.reportConfig = reportConfig;
    }
    return this;
  }

  public ReportConfig getHashConfig() {
    return hashConfig;
  }

  public RunAndReportRequest setHashConfig(ReportConfig hashConfig) {
    if (hashConfig != null) {
      this.hashConfig = hashConfig;
    }
    return this;
  }

  public UUID getRunId() {
    return runId;
  }

  public RunAndReportRequest setRunId(UUID runId) {
    if (runId != null) {
      this.runId = runId;
    }
    return this;
  }
}
