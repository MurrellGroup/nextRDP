interface MetricProps {
  label: string;
  value: string;
  detail?: string;
}

export function Metric({ label, value, detail }: MetricProps) {
  return (
    <div className="metric">
      <span>{label}</span>
      <strong>{value}</strong>
      {detail ? <small>{detail}</small> : null}
    </div>
  );
}
