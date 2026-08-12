import {
  Check,
  FileDown,
  FlaskConical,
  ListChecks,
  ScanSearch,
  SlidersHorizontal,
} from "lucide-react";

import type { WorkflowStep } from "../lib/types";

const steps: Array<{
  id: WorkflowStep;
  label: string;
  caption: string;
  icon: typeof FlaskConical;
}> = [
  { id: "dataset", label: "Dataset", caption: "Load and curate", icon: FlaskConical },
  { id: "settings", label: "Settings", caption: "Analysis design", icon: SlidersHorizontal },
  { id: "scan", label: "Scan", caption: "Detect signals", icon: ScanSearch },
  { id: "review", label: "Review", caption: "Test the hypothesis", icon: ListChecks },
  { id: "export", label: "Export", caption: "Save the analysis", icon: FileDown },
];

interface WorkflowNavProps {
  current: WorkflowStep;
  enabled: Set<WorkflowStep>;
  completed: Set<WorkflowStep>;
  onSelect: (step: WorkflowStep) => void;
}

export function WorkflowNav({ current, enabled, completed, onSelect }: WorkflowNavProps) {
  return (
    <nav className="workflow-nav" aria-label="RDP analysis workflow">
      <ol>
        {steps.map((step, index) => {
          const Icon = step.icon;
          const isCurrent = current === step.id;
          const isEnabled = enabled.has(step.id);
          const isComplete = completed.has(step.id);
          return (
            <li key={step.id}>
              <button
                className={`workflow-step${isCurrent ? " is-current" : ""}${
                  isComplete ? " is-complete" : ""
                }`}
                type="button"
                onClick={() => onSelect(step.id)}
                disabled={!isEnabled}
                aria-current={isCurrent ? "step" : undefined}
              >
                <span className="workflow-index" aria-hidden="true">
                  {isComplete ? <Check size={15} strokeWidth={2.5} /> : index + 1}
                </span>
                <Icon className="workflow-icon" size={18} aria-hidden="true" />
                <span>
                  <strong>{step.label}</strong>
                  <small>{step.caption}</small>
                </span>
              </button>
            </li>
          );
        })}
      </ol>
    </nav>
  );
}
