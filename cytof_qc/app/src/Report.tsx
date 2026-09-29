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

function ClusterProfiles({ report }: { report: Report }) {
  const channels = report.phenotypic_channels;
  const clusters = Object.keys(report.cluster_profiles).sort(
    (a, b) => Number(a) - Number(b)
  );
  // per-channel range normalization so each column reads as its own scale
  const maxOf: Record<string, number> = {};
  for (const ch of channels) {
    maxOf[ch] = Math.max(
      ...clusters.map((c) => report.cluster_profiles[c][ch] ?? 0)
    ) || 1;
  }
  return (
    <section>
      <h3>Cluster marker profiles (median arcsinh)</h3>
      <table className="heatmap">
        <thead>
          <tr>
            <th>cluster</th>
            {channels.map((ch) => (
              <th key={ch}>{ch}</th>
            ))}
          </tr>
        </thead>
        <tbody>
          {clusters.map((c) => (
            <tr key={c}>
              <td>
                {c} ({report.cluster_sizes[c]?.toLocaleString()})
              </td>
              {channels.map((ch) => {
                const v = report.cluster_profiles[c][ch] ?? 0;
                const t = Math.min(1, v / maxOf[ch]);
                return (
                  <td
                    key={ch}
                    style={{
                      background: `rgba(78, 121, 167, ${(0.08 + 0.75 * t).toFixed(2)})`,
                    }}
                  >
                    {v.toFixed(1)}
                  </td>
                );
              })}
            </tr>
          ))}
        </tbody>
      </table>
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
        {report.timing && (
          <p className="scope">
            analyzed in {(report.timing.analyze_ms / 1000).toFixed(1)}s
          </p>
        )}
      </header>

      <section>
        <h3>UMAP by cluster</h3>
        <Scatter
          embedding={report.embedding}
          clusters={report.embedding_clusters}
        />
      </section>

      <ClusterProfiles report={report} />

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
