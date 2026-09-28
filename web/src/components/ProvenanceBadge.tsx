interface Props {
  source: string | undefined
}

export function provenanceLabel(source: string | undefined): string {
  const s = (source ?? '').trim()
  if (!s || s === 'unverified') return 'Unverified provenance'
  if (/synthetic|simulated|fixture/i.test(s)) return `Synthetic data (${s})`
  return `Recorded source: ${s}`
}

export function ProvenanceBadge({ source }: Props) {
  const s = (source ?? '').trim()
  const kind = !s || s === 'unverified'
    ? 'unverified'
    : /synthetic|simulated|fixture/i.test(s)
      ? 'synthetic'
      : 'recorded'
  return (
    <span className={`prov-badge prov-${kind}`} title={`data_source: ${s || 'not recorded'}`}>
      {provenanceLabel(source)}
    </span>
  )
}
