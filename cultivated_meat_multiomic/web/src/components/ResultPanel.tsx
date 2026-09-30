import type { Prediction } from '../api'
import { ProvenanceBadge } from './ProvenanceBadge'

const STATE_COLORS: Record<string, string> = {
  expansion_competent: 'var(--state-green)',
  committed: 'var(--state-amber)',
  terminal: 'var(--state-red)',
}

export function prettyState(s: string): string {
  return s.replace(/_/g, ' ').replace(/\b\w/g, (c) => c.toUpperCase())
}

interface Props {
  result: Prediction | null
  error: string | null
  loading: boolean
}

export function ResultPanel({ result, error, loading }: Props) {
  const entries = result
    ? Object.entries(result.probabilities).sort((a, b) => b[1] - a[1])
    : []

  return (
    <section className="card result-card">
      <div className="card-head">
        <h2>Result</h2>
        {result && <ProvenanceBadge source={result.data_source} />}
      </div>

      {error && <p className="field-error">{error}</p>}

      {!result && !error && (
        <p className="hint">{loading ? 'Predicting…' : 'Run a prediction to see the state distribution.'}</p>
      )}

      {result && (
        <>
          <div
            className="state-banner"
            style={{ borderColor: STATE_COLORS[result.prediction] ?? 'var(--accent)' }}
          >
            <span className="state-label">Predicted state</span>
            <span className="state-name" style={{ color: STATE_COLORS[result.prediction] ?? 'var(--accent)' }}>
              {prettyState(result.prediction)}
            </span>
            <span className="state-conf">confidence {(result.confidence * 100).toFixed(1)}%</span>
          </div>

          <div className="prob-list">
            {entries.map(([state, p]) => (
              <div key={state} className="prob-row">
                <span className="prob-name">{prettyState(state)}</span>
                <div className="prob-track">
                  <div
                    className="prob-fill"
                    style={{
                      width: `${(p * 100).toFixed(1)}%`,
                      background: STATE_COLORS[state] ?? 'var(--accent)',
                    }}
                  />
                </div>
                <span className="prob-value">{(p * 100).toFixed(1)}%</span>
              </div>
            ))}
          </div>
        </>
      )}
    </section>
  )
}
