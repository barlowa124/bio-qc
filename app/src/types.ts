export interface ChannelQc {
  n_events: number;
  negative_fraction: number;
  zero_fraction: number;
  median: number;
  flag_negative: boolean;
  flag_zero: boolean;
}

export interface DriftQc {
  bin_medians: number[];
  median_drift: number;
  flag_drift: boolean;
}

export interface Report {
  n_events: number;
  n_channels: number;
  phenotypic_channels: string[];
  qc_only_channels: string[];
  n_clusters_found: number;
  cluster_sizes: Record<string, number>;
  cluster_profiles: Record<string, Record<string, number>>;
  embedding: number[][];
  embedding_clusters: (string | number)[];
  channels: Record<string, ChannelQc>;
  acquisition_drift: Record<string, DriftQc> | Record<string, Record<string, DriftQc>>;
  timing?: { load_ms: number; analyze_ms: number };
  scope: string;
}
