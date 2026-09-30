import { useCallback, useEffect, useState } from 'react'
import { api, ApiError, type Prediction, type ServiceInfo } from './api'
import { GenePanel } from './components/GenePanel'
import { ResultPanel } from './components/ResultPanel'
import { ProvenanceBadge } from './components/ProvenanceBadge'

type ConnState = 'loading' | 'online' | 'offline'

const POLL_INTERVAL_MS = 30_000

export default function App() {
  const [info, setInfo] = useState<ServiceInfo | null>(null)
  const [conn, setConn] = useState<ConnState>('loading')
  const [values, setValues] = useState<number[]>([])
  const [result, setResult] = useState<Prediction | null>(null)
  const [error, setError] = useState<string | null>(null)
  const [busy, setBusy] = useState(false)

  const probe = useCallback(async () => {
    try {
      const [i] = await Promise.all([api.info(), api.health()])
      setInfo(i)
      setValues((v) => (v.length === i.genes.length ? v : i.genes.map(() => 0)))
      setConn('online')
    } catch {
      setConn('offline')
    }
  }, [])

  useEffect(() => {
    probe()
    const t = setInterval(probe, POLL_INTERVAL_MS)
    return () => clearInterval(t)
  }, [probe])

  function setGene(i: number, v: number) {
    setValues((prev) => {
      const next = [...prev]
      next[i] = v
      return next
    })
    setError(null)
  }

  async function predict() {
    setBusy(true)
    setError(null)
    try {
      setResult(await api.predict(values))
    } catch (e) {
      setResult(null)
      setError(e instanceof ApiError ? e.message : 'Prediction failed.')
      if (e instanceof ApiError && e.status === 0) setConn('offline')
    } finally {
      setBusy(false)
    }
  }

  const offline = conn !== 'online'

  return (
    <div className="app">
      <header className="topbar">
        <div className="brand">
          <span className="eyebrow">Manufacturing QC · 30-Gene Panel</span>
          <h1>Cultivated Meat QC Panel</h1>
          <p className="sub">Manufacturing-readiness state prediction for cultivated muscle tissue</p>
        </div>
        <div className={`conn conn-${conn}`}>
          <span className="conn-dot" />
          {conn === 'loading' ? 'Connecting…' : conn === 'online' ? 'API online' : 'API offline'}
        </div>
      </header>

      {offline && conn === 'offline' && (
        <div className="banner">
          Inference API is not reachable. Start it with <code>cd api && python app.py</code> (binds 127.0.0.1:5000), then wait for the status pill.
        </div>
      )}

      {info && (
        <div className="meta-strip">
          <span><strong>{info.genes.length}</strong> genes</span>
          <span><strong>{info.states.length}</strong> states</span>
          <ProvenanceBadge source={info.data_source} />
        </div>
      )}

      <main className="layout">
        {info && (
          <GenePanel
            genes={info.genes}
            values={values}
            disabled={offline || busy}
            onChange={setGene}
            onBulkSet={(v) => { setValues(v); setError(null) }}
            onReset={() => { setValues(info.genes.map(() => 0)); setResult(null); setError(null) }}
            onPredict={predict}
          />
        )}
        <ResultPanel result={result} error={error} loading={busy} />
      </main>

      <footer className="foot">
        Methods demonstration. Predictions describe melanoma-derived pseudo-bulk clusters,
        not cultivated-meat cell states.
      </footer>
    </div>
  )
}
