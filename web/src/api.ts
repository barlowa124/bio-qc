export interface ServiceInfo {
  service: string
  genes: string[]
  states: string[]
  data_source: string
}

export interface Health {
  status: string
  model_accuracy: number | null
  data_source: string
}

export interface Prediction {
  prediction: string
  probabilities: Record<string, number>
  confidence: number
  data_source: string
}

export class ApiError extends Error {
  status: number
  constructor(status: number, message: string) {
    super(message)
    this.status = status
    this.name = 'ApiError'
  }
}

const BASE = (import.meta.env.VITE_API_BASE ?? '/api').replace(/\/$/, '')

async function request<T>(path: string, init?: RequestInit): Promise<T> {
  let res: Response
  try {
    res = await fetch(`${BASE}${path}`, init)
  } catch {
    throw new ApiError(0, 'API unreachable')
  }
  const body = await res.json().catch(() => ({}))
  if (!res.ok) {
    throw new ApiError(res.status, body.error ?? `HTTP ${res.status}`)
  }
  return body as T
}

export const api = {
  info: () => request<ServiceInfo>('/'),
  health: () => request<Health>('/health'),
  predict: (expression: number[]) =>
    request<Prediction>('/predict', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ expression }),
    }),
}
