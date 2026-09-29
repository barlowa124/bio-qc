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
      // analysis runs as a job — poll until the report is ready
      const jobId: string = body.job_id;
      for (let i = 0; i < 150; i++) {
        await new Promise((r) => setTimeout(r, 2000));
        const jr = await fetch(`/api/jobs/${jobId}`);
        const jb = await jr.json();
        if (jb.status === 'running') continue;
        if (jb.status === 'error') throw new Error(jb.error);
        setReport(jb);
        return;
      }
      throw new Error('analysis timed out after 5 minutes');
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
        {busy ? 'analysis running — large uploads take a minute or two' : 'choose files'}
      </label>
      <p>
        or{' '}
        <button
          type="button"
          className="link"
          onClick={async () => {
            setError(null);
            setBusy(true);
            try {
              const r = await fetch('/example_report.json');
              if (!r.ok) throw new Error(`HTTP ${r.status}`);
              setReport(await r.json());
            } catch (e) {
              setError(e instanceof Error ? e.message : String(e));
            } finally {
              setBusy(false);
            }
          }}
        >
          load an example report
        </button>{' '}
        (9,222 Levine bone-marrow events, precomputed)
      </p>
      {error && <p className="err">{error}</p>}
      {report && <ReportView report={report} />}
    </main>
  );
}
