# TCRtoolkit (Integrated) — Design

Merged pipeline: bulk TCR analysis + single-cell TCR analysis sharing one bulk engine.
Full spec in `IMPLEMENTATION_SPEC.md` (same folder). This file is the visual design.

Legend: solid = active path · dashed = skipped when GEX object absent · ⭐ = shared bulk engine
(unchanged) · 🧬 = single-cell-specific · 🔀 = bridge (schema conversion).

---

## 1. Top-level modality dispatch

```mermaid
flowchart TD
    RUN(["nextflow run . --mode ?"]) --> D{"--mode<br/>(default: bulk)"}
    D -->|bulk| BULK["BULK route<br/>airr / adaptive input"]
    D -->|singlecell| SC["SINGLE-CELL route<br/>cellranger + sample sheet<br/>(+ optional GEX object)"]

    BULK --> ENG
    SC --> ENG["⭐ Shared tcrtoolkit bulk engine<br/>ANNOTATE · SAMPLE · PATIENT · COMPARE"]

    style ENG fill:#1f6f43,stroke:#0d3,color:#fff
```

---

## 2. Full integrated pipeline (the spine)

```mermaid
flowchart TD
    %% ---------- inputs ----------
    IN["🧬 Cell Ranger VDJ contigs<br/>+ sample sheet<br/>(+ optional GEX annotated object)"]

    IN --> VQ["🧬 VDJ_QC<br/>contig QC + repertoire QC tables"]
    VQ --> GEX{"GEX annotated<br/>object provided?"}

    %% ---------- head: GEX-only ----------
    GEX -->|yes| TI["🧬 TCELL_INTEGRATION<br/>filter T-cells · per-cell object"]
    TI --> PBp["🔀 PSEUDOBULK (by phenotype)<br/>→ canonical AIRR schema"]
    GEX -->|no| PBs["🔀 PSEUDOBULK (single group)<br/>→ canonical AIRR schema"]

    PBp --> QC
    PBs --> QC["🧬 PSEUDOBULK_QC<br/>gate n_clones / n_cells"]

    %% ---------- shared bulk engine ----------
    QC --> AFC["⭐ ANNOTATE_FROM_CONCAT<br/>sort · dedup · OLGA"]
    AFC --> SA["⭐ SAMPLE<br/>Shannon/Simpson/convergence · tcrdist3"]
    AFC --> PA["⭐ PATIENT<br/>GIANA + GLIPH2 (emits clusters)"]
    AFC --> CO["⭐ COMPARE"]

    %% ---------- tail: GEX-only ----------
    PA --> CTS["🔀 CLUSTER_TO_SC<br/>map clusters → cells → enriched Seurat"]
    SA --> CTS
    CTS -.GEX only.-> CG["🧬 CONGA<br/>(needs gene expression)"]
    CTS --> RE

    %% repertoire runs in both, adaptively
    SA --> RE["🧬 REPERTOIRE (adaptive)<br/>cell-level if GEX · clonotype-level if not"]
    CG --> MS
    RE --> MS["🧬 MASTER_SUMMARY<br/>conditional sections"]
    CO --> MS

    %% ---------- styling ----------
    style AFC fill:#1f6f43,stroke:#0d3,color:#fff
    style SA fill:#1f6f43,stroke:#0d3,color:#fff
    style PA fill:#1f6f43,stroke:#0d3,color:#fff
    style CO fill:#1f6f43,stroke:#0d3,color:#fff
    style CG stroke-dasharray: 5 5
```

---

## 3. Two routes overlaid (what runs vs. what is skipped)

```mermaid
flowchart LR
    subgraph SHARED ["Runs in BOTH routes (identical code)"]
      direction TB
      A1[VDJ_QC] --> A2[PSEUDOBULK] --> A3[PSEUDOBULK_QC]
      A3 --> A4["⭐ ANNOTATE_FROM_CONCAT"] --> A5["⭐ SAMPLE"]
      A4 --> A6["⭐ PATIENT: GIANA/GLIPH2"] --> A7["⭐ COMPARE"]
      A5 --> A8["REPERTOIRE (adaptive)"]
    end

    subgraph GEXONLY ["Adds ONLY when GEX object present"]
      direction TB
      B1[TCELL_INTEGRATION] -.-> B2["PSEUDOBULK: by phenotype"]
      B3[CLUSTER_TO_SC] -.-> B4[CONGA]
      B5["REPERTOIRE: cell-level + UMAP"]
      B6["MASTER_SUMMARY: full"]
    end

    SHARED --> MS2["MASTER_SUMMARY<br/>(full if GEX, reduced if not)"]
    GEXONLY --> MS2
```

---

## 4. Where the code lives (module → source → container)

```mermaid
flowchart TD
    subgraph SCRATCH ["🧬 single-cell modules — container: scratch-tcr"]
      VDJ_QC2[VDJ_QC]
      TI2[TCELL_INTEGRATION]
      CONGA2[CONGA]
      CONS2[CONSENSUS_CLUSTERING]
      REP2[REPERTOIRE]
      MAS2[MASTER_SUMMARY]
    end

    subgraph BRIDGES ["🔀 bridges — schema conversion"]
      SCC2[SC_TO_CDR3 = pseudobulk]
      VTB2[VDJ_TO_BULK]
      CTS2[CLUSTER_TO_SC]
    end

    subgraph ENGINE ["⭐ shared bulk engine — container: tcrtoolkit:main (UNCHANGED)"]
      ANN2["ANNOTATE / ANNOTATE_FROM_CONCAT*"]
      SAM2[SAMPLE]
      PAT2["PATIENT*"]
      COM2[COMPARE]
      GIA2[GIANA_CALC]
      GLI2[GLIPH2_TURBOGLIPH]
      TCR2[TCRDIST3_MATRIX]
    end

    BRIDGES --> ENGINE
    SCRATCH --> BRIDGES

    note["* = additive edits only:<br/>+ ANNOTATE_FROM_CONCAT entrypoint<br/>+ PATIENT cluster emits<br/>+ optional patient_col"]
    ENGINE -.-> note
    style note fill:#332,stroke:#a80,color:#fff
```

---

## 5. Bulk-impact summary

```mermaid
flowchart LR
    subgraph FROZEN ["Bulk: unchanged behavior"]
      F1[airr input] --- F2[adaptive input]
      F3[ANNOTATE] --- F4[SAMPLE] --- F5[COMPARE]
      F6[GIANA / GLIPH2 / tcrdist modules]
    end

    subgraph ADDED ["Bulk: additive edits (no behavior change)"]
      G1["+ ANNOTATE_FROM_CONCAT (new entrypoint)"]
      G2["+ PATIENT emits + optional patient_col"]
    end

    subgraph MOVED ["Relocated OUT of bulk → into single-cell"]
      H1["cellranger pseudobulk"]
      H2["phenotype-stratified pseudobulk"]
      H3["sobject_gex param"]
    end
```
