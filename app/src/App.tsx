import { useState } from 'react';
import { Report } from './types';
import { ReportView } from './Report';

export function App() {
  const [report, setReport] = useState<Report | null>(null);
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState<string | null>(null);

  async function onFiles(files: FileList | null) {
    if (!files || files.length === 0) return;
    setBusy(true);
    setError(null);
    try {
      const fd = new FormData();
      for (const f of Array.from(files)) fd.append('files', f);
      const resp = await fetch('/api/analyze', { method: 'POST', body: fd });
      const body = await resp.json();
      if (!resp.ok) throw new Error(body.error || `HTTP ${resp.status}`);
      setReport(body);
    } catch (e) {
      setError(e instanceof Error ? e.message : String(e));
    } finally {
      setBusy(false);
    }
  }

  return (
    <main>
      <h1>cytof-qc</h1>
      <p className="tag">
        Upload CyTOF events (.fcs or .csv) for channel QC, acquisition-drift
        flags, and Leiden clustering.
      </p>
      <label className="drop">
        <input
          type="file"
          multiple
          accept=".fcs,.csv"
          onChange={(e) => onFiles(e.target.files)}
        />
        {busy ? 'analyzing…' : 'choose files'}
      </label>
      {error && <p className="err">{error}</p>}
      {report && <ReportView report={report} />}
    </main>
  );
}
