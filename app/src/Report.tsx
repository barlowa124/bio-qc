import { Report, DriftQc } from './types';
import { Scatter } from './Scatter';

function DriftTable({ drift }: { drift: Report['acquisition_drift'] }) {
  // drift may be flat {channel: DriftQc} or nested {batch: {channel: DriftQc}}
  const flat: Record<string, DriftQc & { batch?: string }> = {};
  for (const [k, v] of Object.entries(drift)) {
    if (v && typeof v === 'object' && 'median_drift' in v) {
      flat[k] = v as DriftQc;
    } else {
      for (const [ch, inner] of Object.entries(v as Record<string, DriftQc>)) {
        flat[`${k} / ${ch}`] = { ...inner, batch: k };
      }
    }
  }
  const flagged = Object.entries(flat).filter(([, d]) => d.flag_drift);
  return (
    <section>
      <h3>Acquisition drift</h3>
      {flagged.length === 0 ? (
        <p className="ok">No channel exceeds the drift gate.</p>
      ) : (
        <table>
          <thead>
            <tr><th>channel</th><th>median drift (arcsinh)</th></tr>
          </thead>
          <tbody>
            {flagged.map(([name, d]) => (
              <tr key={name} className="flag">
                <td>{name}</td>
                <td>{d.median_drift.toFixed(3)}</td>
              </tr>
            ))}
          </tbody>
        </table>
      )}
    </section>
  );
}

export function ReportView({ report }: { report: Report }) {
  const channelRows = Object.entries(report.channels);
  const flaggedChannels = channelRows.filter(
    ([, c]) => c.flag_negative || c.flag_zero
  );
  return (
    <div className="report">
      <header>
        <h2>
          {report.n_events.toLocaleString()} events ·{' '}
          {report.n_clusters_found} clusters · {report.n_channels} markers
        </h2>
        <p className="scope">{report.scope}</p>
      </header>

      <section>
        <h3>UMAP by cluster</h3>
        <Scatter
          embedding={report.embedding}
          clusters={report.embedding_clusters}
        />
      </section>

      <DriftTable drift={report.acquisition_drift} />

      <section>
        <h3>Channel QC</h3>
        <table>
          <thead>
            <tr>
              <th>channel</th><th>neg %</th><th>zero %</th><th>median</th>
            </tr>
          </thead>
          <tbody>
            {channelRows.map(([name, c]) => (
              <tr
                key={name}
                className={c.flag_negative || c.flag_zero ? 'flag' : ''}
              >
                <td>{name}</td>
                <td>{(c.negative_fraction * 100).toFixed(1)}</td>
                <td>{(c.zero_fraction * 100).toFixed(1)}</td>
                <td>{c.median.toFixed(2)}</td>
              </tr>
            ))}
          </tbody>
        </table>
        {flaggedChannels.length === 0 && (
          <p className="ok">All channels inside QC tolerances.</p>
        )}
      </section>
    </div>
  );
}
