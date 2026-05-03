package neqsim.process.processmodel;

import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;
import org.junit.jupiter.api.Test;

class ProcessJsonValidatorTest {
  @Test
  void validJsonShouldPass() {
    String json = "{\"process\":[{\"type\":\"Stream\",\"name\":\"feed\"},{\"type\":\"Separator\",\"name\":\"sep\",\"inlet\":\"feed\"}]}";
    ProcessJsonValidator.ValidationReport report = ProcessJsonValidator.validate(json);
    assertTrue(report.isValid());
    assertTrue(report.getWarningCount() == 0);
  }

  @Test
  void missingProcessShouldFail() {
    ProcessJsonValidator.ValidationReport report = ProcessJsonValidator.validate("{}");
    assertFalse(report.isValid());
  }

  @Test
  void duplicateNamesShouldFail() {
    String json = "{\"process\":[{\"type\":\"Stream\",\"name\":\"A\"},{\"type\":\"Stream\",\"name\":\"A\"}]}";
    ProcessJsonValidator.ValidationReport report = ProcessJsonValidator.validate(json);
    assertFalse(report.isValid());
  }

  @Test
  void unknownStreamReferenceShouldWarnOnly() {
    String json = "{\"process\":[{\"type\":\"Separator\",\"name\":\"sep\",\"inlet\":\"missing\"}]}";
    ProcessJsonValidator.ValidationReport report = ProcessJsonValidator.validate(json);
    assertTrue(report.isValid());
    assertTrue(report.getWarningCount() > 0);
  }

  @Test
  void processSystemValidateJsonMethodShouldWork() {
    String json = "{\"process\":[{\"type\":\"Stream\",\"name\":\"feed\"}]}";
    ProcessJsonValidator.ValidationReport report = ProcessSystem.validateJson(json);
    assertTrue(report.isValid());
  }

  @Test
  void nullJsonShouldFail() {
    ProcessJsonValidator.ValidationReport report = ProcessJsonValidator.validate(null);
    assertFalse(report.isValid());
  }

  @Test
  void malformedJsonShouldFail() {
    ProcessJsonValidator.ValidationReport report = ProcessJsonValidator.validate("{bad");
    assertFalse(report.isValid());
  }

  @Test
  void missingTypeShouldFail() {
    String json = "{\"process\":[{\"name\":\"feed\"}]}";
    ProcessJsonValidator.ValidationReport report = ProcessJsonValidator.validate(json);
    assertFalse(report.isValid());
  }

  @Test
  void missingNameShouldFail() {
    String json = "{\"process\":[{\"type\":\"Stream\"}]}";
    ProcessJsonValidator.ValidationReport report = ProcessJsonValidator.validate(json);
    assertFalse(report.isValid());
  }

  @Test
  void unknownReferenceInStreamsObjectShouldWarn() {
    String json = "{\"process\":[{\"type\":\"Separator\",\"name\":\"sep\",\"streams\":{\"inlet\":\"ghost\"}}]}";
    ProcessJsonValidator.ValidationReport report = ProcessJsonValidator.validate(json);
    assertTrue(report.isValid());
    assertTrue(report.getWarningCount() > 0);
  }

}
