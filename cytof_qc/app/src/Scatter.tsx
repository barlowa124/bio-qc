import { useEffect, useRef } from 'react';

// Categorical palette, colorblind-checked up to ~10 clusters; beyond that
// it cycles, which is acceptable for a QC view where identity matters
// less than separation.
const PALETTE = [
  '#4e79a7', '#f28e2b', '#e15759', '#76b7b2', '#59a14f',
  '#edc948', '#b07aa1', '#ff9da7', '#9c755f', '#bab0ac',
];

export function Scatter({
  embedding,
  clusters,
}: {
  embedding: number[][];
  clusters: (string | number)[];
}) {
  const ref = useRef<HTMLCanvasElement>(null);

  useEffect(() => {
    const cv = ref.current;
    if (!cv || embedding.length === 0) return;
    const dpr = window.devicePixelRatio || 1;
    const w = cv.clientWidth;
    const h = cv.clientHeight;
    cv.width = w * dpr;
    cv.height = h * dpr;
    const g = cv.getContext('2d')!;
    g.scale(dpr, dpr);
    g.fillStyle = '#0b0f14';
    g.fillRect(0, 0, w, h);

    let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
    for (const [x, y] of embedding) {
      if (x < minX) minX = x;
      if (x > maxX) maxX = x;
      if (y < minY) minY = y;
      if (y > maxY) maxY = y;
    }
    const pad = 12;
    const sx = (w - 2 * pad) / ((maxX - minX) || 1);
    const sy = (h - 2 * pad) / ((maxY - minY) || 1);
    const s = Math.min(sx, sy);

    const keys = [...new Set(clusters.map(String))];
    const colorOf = (c: string | number) =>
      PALETTE[keys.indexOf(String(c)) % PALETTE.length];

    const r = embedding.length > 5000 ? 1.2 : 2;
    for (let i = 0; i < embedding.length; i++) {
      const [x, y] = embedding[i];
      g.fillStyle = colorOf(clusters[i]);
      g.beginPath();
      g.arc(pad + (x - minX) * s, pad + (y - minY) * s, r, 0, Math.PI * 2);
      g.fill();
    }
  }, [embedding, clusters]);

  return <canvas ref={ref} className="scatter" />;
}
