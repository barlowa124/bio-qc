from html import escape


def data_source(artifact):
    value = artifact.get("data_source") if isinstance(artifact, dict) else None
    return value.strip() if isinstance(value, str) and value.strip() else "unverified"


def is_synthetic(artifact):
    source = data_source(artifact).lower()
    return any(word in source for word in ("synthetic", "simulated", "fixture"))


def provenance_label(artifact):
    source = data_source(artifact)
    if is_synthetic(artifact):
        return f"Synthetic data ({source})"
    if source == "unverified":
        return "Unverified provenance (source not recorded)"
    return f"Recorded source: {source}"


def annotate_figure(fig, artifact):
    label = provenance_label(artifact)
    title = fig.layout.title.text or ""
    suffix = f"<br><sup>{escape(label)}</sup>"
    if not title.endswith(suffix):
        fig.update_layout(title_text=title + suffix)
    meta = dict(fig.layout.meta) if isinstance(fig.layout.meta, dict) else {}
    meta.update(data_source=data_source(artifact), synthetic=is_synthetic(artifact))
    fig.update_layout(meta=meta)
    return fig
