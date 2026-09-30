import { useRef, useState } from 'react'

interface Props {
  genes: string[]
  values: number[]
  disabled: boolean
  onChange: (index: number, value: number) => void
  onBulkSet: (values: number[]) => void
  onReset: () => void
  onPredict: () => void
}

export function parseExpression(text: string, expected: number): number[] | string {
  const tokens = text.split(/[\s,;\t]+/).filter(Boolean)
  const numeric = tokens.filter((t) => !Number.isNaN(Number(t)))
  if (numeric.length === 0) return 'No numeric values found.'
  if (numeric.length !== expected)
    return `Expected ${expected} values, found ${numeric.length}.`
  return numeric.map(Number)
}

export function GenePanel({ genes, values, disabled, onChange, onBulkSet, onReset, onPredict }: Props) {
  const fileRef = useRef<HTMLInputElement>(null)
  const [pasteOpen, setPasteOpen] = useState(false)
  const [pasteText, setPasteText] = useState('')
  const [bulkError, setBulkError] = useState<string | null>(null)

  function applyBulk(text: string) {
    const parsed = parseExpression(text, genes.length)
    if (typeof parsed === 'string') {
      setBulkError(parsed)
    } else {
      setBulkError(null)
      onBulkSet(parsed)
      setPasteText('')
      setPasteOpen(false)
    }
  }

  function onFile(e: React.ChangeEvent<HTMLInputElement>) {
    const f = e.target.files?.[0]
    if (!f) return
    f.text().then((text) => {
      const lines = text.split(/\r?\n/).filter((l) => l.trim())
      const rows = lines.filter((l) => l.split(/[\s,;\t]+/).some((t) => !Number.isNaN(Number(t)) && t !== ''))
      if (rows.length === 0) setBulkError('No numeric row found in CSV.')
      else applyBulk(rows[0])
    })
    e.target.value = ''
  }

  return (
    <section className="card">
      <div className="card-head">
        <h2>Expression input</h2>
        <div className="card-actions">
          <button type="button" className="btn-ghost" onClick={() => setPasteOpen((v) => !v)}>
            Paste values
          </button>
          <button type="button" className="btn-ghost" onClick={() => fileRef.current?.click()}>
            Upload CSV row
          </button>
          <input ref={fileRef} type="file" accept=".csv,.txt" hidden onChange={onFile} />
        </div>
      </div>
      <p className="hint">Log TPM values for the {genes.length}-gene panel.</p>

      {pasteOpen && (
        <div className="paste-block">
          <textarea
            value={pasteText}
            onChange={(e) => setPasteText(e.target.value)}
            placeholder={`${genes.length} comma- or space-separated values, in panel order`}
            rows={3}
          />
          <div className="paste-actions">
            <button type="button" className="btn-ghost" onClick={() => applyBulk(pasteText)}>
              Fill inputs
            </button>
            {bulkError && <span className="field-error">{bulkError}</span>}
          </div>
        </div>
      )}
      {bulkError && !pasteOpen && <p className="field-error">{bulkError}</p>}

      <div className="gene-grid">
        {genes.map((g, i) => (
          <label key={g} className="gene-field">
            <span className="gene-name">{g}</span>
            <input
              type="number"
              step="0.1"
              value={Number.isFinite(values[i]) ? values[i] : 0}
              disabled={disabled}
              onChange={(e) => onChange(i, Number(e.target.value))}
            />
          </label>
        ))}
      </div>

      <div className="card-foot">
        <button type="button" className="btn-primary" onClick={onPredict} disabled={disabled}>
          Predict state
        </button>
        <button type="button" className="btn-ghost" onClick={onReset} disabled={disabled}>
          Reset
        </button>
      </div>
    </section>
  )
}
