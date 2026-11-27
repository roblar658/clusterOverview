// =============================================================================
// IEEE Transactions Two-Column Mathematical Monograph Template (Typst 0.15)
// Created by Gemini, AI made by Google
// =============================================================================

#let ieee_navy = rgb("#0f294a")
#let ieee_blue = rgb("#1d4ed8")
#let ieee_slate = rgb("#475569")
#let ieee_muted = rgb("#64748b")
#let ieee_border = rgb("#cbd5e1")
#let ieee_bg_orig = rgb("#f8fafc")
#let ieee_bg_cite = rgb("#f8f9fa")
#let ieee_bg_alg = rgb("#fbfcfe")

#let thm_counter = counter("ieee_thm")
#let lem_counter = counter("ieee_lem")
#let prop_counter = counter("ieee_prop")
#let cor_counter = counter("ieee_cor")
#let def_counter = counter("ieee_def")
#let alg_counter = counter("ieee_alg")

#let ieee_doc(
  title: "Mathematical Foundations, Proofs, and Applications of 200 Clustering Algorithms",
  subtitle: "An IEEE-Formatted Reference Monograph with Original Proofs, Verified Classical Derivations, and Cross-Paradigm Comparisons",
  author: "Gemini, AI made by Google",
  affiliation: "Created by Gemini, AI made by Google  •  ClusterOverview Project (2026)",
  doc,
) = {
  set document(
    title: title,
    author: author,
  )

  set page(
    paper: "us-letter",
    margin: (top: 0.72in, bottom: 0.82in, left: 0.60in, right: 0.60in),
    columns: 2,
    header-ascent: 24%,
    footer-descent: 28%,
    header: context {
      let pg = counter(page).get().first()
      if pg > 1 {
        set text(size: 7.5pt, fill: ieee_slate, font: ("Times New Roman", "Libertinus Serif"))
        let is_even = calc.rem(pg, 2) == 0
        if is_even {
          grid(
            columns: (auto, 1fr),
            align: (left, right),
            [#pg],
            [IEEE TRANSACTIONS MONOGRAPH — 200 CLUSTERING ALGORITHMS (CREATED BY GEMINI, AI MADE BY GOOGLE)],
          )
        } else {
          let before_h = query(selector(heading.where(level: 2)).before(here()))
          let cur = if before_h.len() > 0 { before_h.last().body } else { [MATHEMATICAL FOUNDATIONS OF CLUSTERING] }
          grid(
            columns: (1fr, auto),
            align: (left, right),
            [#upper(cur)],
            [#pg],
          )
        }
        v(2pt)
        line(length: 100%, stroke: 0.45pt + ieee_border)
      }
    },
    footer: context {
      let pg = counter(page).get().first()
      if pg > 1 {
        set text(size: 7pt, fill: ieee_muted, font: ("Times New Roman", "Libertinus Serif"))
        line(length: 100%, stroke: 0.35pt + ieee_border)
        v(2pt)
        grid(
          columns: (1fr, auto, 1fr),
          align: (left, center, right),
          [Created by Gemini, AI made by Google],
          [ClusterOverview Mathematical Compendium],
          [Page #pg],
        )
      }
    },
  )

  set columns(gutter: 0.24in)

  set text(
    font: ("Times New Roman", "Libertinus Serif", "New Computer Modern"),
    size: 9.5pt,
    fill: rgb("#0f172a"),
    lang: "en",
    hyphenate: true,
  )

  show math.equation: set text(
    font: ("New Computer Modern Math", "Cambria Math"),
    size: 9.1pt,
  )
  show math.equation.where(block: true): set text(size: 8.7pt)
  show math.equation.where(block: true): set block(above: 0.55em, below: 0.55em, breakable: false)

  set par(
    justify: true,
    first-line-indent: (amount: 1em, all: false),
    leading: 0.53em,
    spacing: 0.56em,
  )

  set list(indent: 0.9em, body-indent: 0.45em, spacing: 0.35em)
  set enum(indent: 0.9em, body-indent: 0.45em, spacing: 0.35em)

  // IEEE Booktabs Table Defaults
  set table(
    stroke: (_, y) => (
      top: if y == 0 { 0.9pt + ieee_navy } else if y == 1 { 0.55pt + ieee_navy } else { 0pt },
      bottom: 0.4pt + ieee_border,
    ),
    inset: (x: 4.5pt, y: 3.5pt),
  )
  show table: set text(size: 7.8pt)
  show table.cell.where(y: 0): set text(weight: "bold", fill: ieee_navy)

  // Heading hierarchy in authentic IEEE style
  show heading.where(level: 1): it => {
    set align(center)
    set par(first-line-indent: 0pt, justify: false)
    text(size: 13.5pt, weight: "bold", tracking: 0.5pt, fill: white, upper(it.body))
  }

  show heading.where(level: 2): it => {
    v(0.75em, weak: true)
    set par(first-line-indent: 0pt)
    set align(center)
    text(size: 10pt, weight: "bold", fill: ieee_navy, upper(it.body))
    v(0.35em, weak: true)
  }

  show heading.where(level: 3): it => {
    v(0.55em, weak: true)
    set par(first-line-indent: 0pt)
    text(size: 9.5pt, weight: "bold", style: "italic", fill: ieee_navy, it.body)
    v(0.22em, weak: true)
  }

  show heading.where(level: 4): it => {
    v(0.35em, weak: true)
    set par(first-line-indent: 0pt)
    text(size: 9pt, weight: "bold", fill: ieee_slate, it.body)
    v(0.15em, weak: true)
  }

  set bibliography(style: "ieee")

  doc
}

// =============================================================================
// Full-Width IEEE Part Opener Page (for Parts I–IX, never a blank stub page)
// =============================================================================
#let part_opener(
  roman: "I",
  title: "Centroid, Medoid, and Partitioning Algorithms",
  alg_range: "Algorithms 1–30",
  summary: [],
  core_principle: [],
  alg_list: [],
  list_heading: none,
) = {
  pagebreak(weak: true)
  place(
    top + center,
    scope: "parent",
    float: true,
    block(
      width: 100%,
      above: 6pt,
      below: 12pt,
      {
        block(
          width: 100%,
          inset: (x: 16pt, y: 16pt),
          fill: ieee_navy,
          radius: 3pt,
          {
            set align(center)
            set par(first-line-indent: 0pt, justify: false)
            text(size: 9pt, weight: "bold", tracking: 1.5pt, fill: rgb("#93c5fd"))[
              PART #roman  •  #upper(alg_range)
            ]
            v(5pt)
            heading(level: 1, outlined: true, numbering: none)[Part #roman: #title]
          },
        )
        v(10pt)
        grid(
          columns: (1fr, 1fr),
          column-gutter: 14pt,
          block(
            width: 100%,
            inset: (left: 9pt, right: 8pt, top: 7pt, bottom: 7pt),
            stroke: (left: 2.2pt + ieee_blue),
            fill: rgb("#f8fafc"),
            {
              set par(first-line-indent: 0pt, justify: true)
              set text(size: 9pt)
              text(weight: "bold", fill: ieee_navy)[Paradigm Overview: ]
              summary
            },
          ),
          block(
            width: 100%,
            inset: (left: 9pt, right: 8pt, top: 7pt, bottom: 7pt),
            stroke: (left: 2.2pt + ieee_slate),
            fill: rgb("#f8f9fa"),
            {
              set par(first-line-indent: 0pt, justify: true)
              set text(size: 9pt)
              text(weight: "bold", fill: ieee_navy)[Unifying Mathematical Principle: ]
              core_principle
            },
          ),
        )
        if alg_list != [] {
          v(12pt)
          block(
            width: 100%,
            stroke: (top: 0.9pt + ieee_navy, bottom: 0.9pt + ieee_navy),
            inset: (x: 8pt, y: 8pt),
            fill: white,
            {
              set par(first-line-indent: 0pt, justify: false, leading: 0.55em)
              text(size: 8.8pt, weight: "bold", fill: ieee_navy)[#(if list_heading != none { list_heading } else [ALGORITHMS COVERED IN PART #roman])]
              v(5pt)
              set text(size: 8.3pt)
              alg_list
            },
          )
        }
      },
    ),
  )
}

// =============================================================================
// Full-width IEEE Algorithm Chapter Banner (starts each algorithm on a clean page)
// =============================================================================
#let alg_entry(
  id: 1,
  name: "Algorithm Name",
  category: "Category I: Centroid, Medoid & Partitioning",
  year: "1982",
  paradigm: "Prototype Minimization",
  complexity: $O(n k d T)$,
  provenance: "Original Derivation + Cited Classical Bounds",
) = {
  pagebreak(weak: true)
  thm_counter.update(0)
  lem_counter.update(0)
  prop_counter.update(0)
  cor_counter.update(0)
  def_counter.update(0)
  alg_counter.update(0)

  place(
    top,
    scope: "parent",
    float: true,
    block(
      width: 100%,
      above: 0pt,
      below: 10pt,
      stroke: (top: 1.8pt + ieee_navy, bottom: 0.9pt + ieee_navy),
      inset: (x: 6pt, y: 7pt),
      fill: rgb("#f8fafc"),
      {
        set par(first-line-indent: 0pt, justify: false)
        grid(
          columns: (auto, 1fr, auto),
          column-gutter: 8pt,
          align: (left + horizon, left + horizon, right + horizon),
          box(
            fill: ieee_navy,
            inset: (x: 6pt, y: 3pt),
            radius: 2pt,
            text(size: 9.5pt, weight: "bold", fill: white)[ALG. \##id],
          ),
          [
            #heading(level: 2, outlined: true, numbering: none)[Algorithm #id: #name]
          ],
          text(size: 8.5pt, style: "italic", fill: ieee_slate)[Introduced: #year],
        )
        v(3pt)
        line(length: 100%, stroke: 0.35pt + ieee_border)
        v(3pt)
        set text(size: 8pt, fill: ieee_slate)
        grid(
          columns: (1.25fr, 1fr, 0.85fr, 1.2fr),
          column-gutter: 6pt,
          [*Family:* #category],
          [*Paradigm:* #paradigm],
          [*Complexity:* #complexity],
          [*Proof Basis:* #provenance],
        )
      },
    ),
  )
}

// =============================================================================
// Full-width IEEE Comparative Chapter Banner (Part IX)
// =============================================================================
#let cmp_entry(
  roman: "I",
  title: "Comparative Chapter Title",
  scope_label: "Cross-Paradigm Theoretical Analysis",
  provenance: "Original Unifying Proofs + Cited Surveys",
) = {
  pagebreak(weak: true)
  thm_counter.update(0)
  lem_counter.update(0)
  prop_counter.update(0)
  cor_counter.update(0)
  def_counter.update(0)
  alg_counter.update(0)

  place(
    top,
    scope: "parent",
    float: true,
    block(
      width: 100%,
      above: 0pt,
      below: 10pt,
      stroke: (top: 1.8pt + ieee_blue, bottom: 0.9pt + ieee_navy),
      inset: (x: 8pt, y: 8pt),
      fill: rgb("#eff6ff"),
      {
        set par(first-line-indent: 0pt, justify: false)
        grid(
          columns: (auto, 1fr),
          column-gutter: 8pt,
          align: (left + horizon, left + horizon),
          box(
            fill: ieee_blue,
            inset: (x: 6pt, y: 3pt),
            radius: 2pt,
            text(size: 9.2pt, weight: "bold", fill: white)[COMP. CH. #roman],
          ),
          [
            #heading(level: 2, outlined: true, numbering: none)[Comparison Chapter #roman: #title]
          ],
        )
        v(3pt)
        line(length: 100%, stroke: 0.35pt + ieee_border)
        v(3pt)
        set text(size: 8pt, fill: ieee_slate)
        grid(
          columns: (1fr, 1fr),
          column-gutter: 8pt,
          [*Scope:* #scope_label],
          [*Proof Basis:* #provenance],
        )
      },
    ),
  )
}

// =============================================================================
// Section Subheaders inside an Algorithm Entry (IEEE Roman & Letter conventions)
// =============================================================================
#let sec(roman, title) = {
  v(0.65em, weak: true)
  set par(first-line-indent: 0pt)
  block(
    sticky: true,
    below: 0.32em,
    {
      set align(center)
      text(size: 9.5pt, weight: "bold", fill: ieee_navy)[#roman. #upper(title)]
    },
  )
}

#let subsec(letter, title) = {
  v(0.45em, weak: true)
  set par(first-line-indent: 0pt)
  block(
    sticky: true,
    below: 0.2em,
    text(size: 9.3pt, weight: "bold", style: "italic", fill: ieee_navy)[#letter. #title],
  )
}

// =============================================================================
// Formal Mathematical Environments with Provenance Tags
// =============================================================================
#let _origin_badge(origin) = {
  let is_orig = lower(origin).contains("original") or lower(origin).contains("constructive") or lower(origin).contains("unified")
  let col = if is_orig { ieee_blue } else { ieee_slate }
  let bg = if is_orig { rgb("#eff6ff") } else { rgb("#f1f5f9") }
  box(
    fill: bg,
    stroke: 0.45pt + col,
    radius: 2pt,
    inset: (x: 3.5pt, y: 1.2pt),
    baseline: 12%,
    text(size: 6.8pt, weight: "bold", fill: col)[#upper(origin)],
  )
}

#let definition(title: none, body) = {
  def_counter.step()
  block(
    width: 100%,
    above: 0.55em,
    below: 0.55em,
    inset: (left: 7pt, right: 5pt, top: 4.5pt, bottom: 4.5pt),
    stroke: (left: 1.8pt + ieee_slate),
    fill: rgb("#f8fafc"),
    breakable: true,
    {
      set par(first-line-indent: 0pt)
      context {
        let n = def_counter.get().first()
        text(weight: "bold", fill: ieee_navy)[Definition #n]
        if title != none { text(weight: "bold", fill: ieee_navy)[ (#title)] }
        text(weight: "bold", fill: ieee_navy)[: ]
      }
      body
    },
  )
}

#let _stmt_box(kind, ctr, title, origin, body) = {
  ctr.step()
  let is_orig = origin != none and (lower(origin).contains("original") or lower(origin).contains("constructive") or lower(origin).contains("unified"))
  let accent = if is_orig { ieee_blue } else { ieee_navy }
  let bg = if is_orig { ieee_bg_orig } else { ieee_bg_cite }
  block(
    width: 100%,
    above: 0.55em,
    below: 0.35em,
    inset: (left: 7pt, right: 5pt, top: 5pt, bottom: 5pt),
    stroke: (left: 2pt + accent, top: 0.35pt + ieee_border, right: 0.35pt + ieee_border, bottom: 0.35pt + ieee_border),
    fill: bg,
    breakable: true,
    {
      set par(first-line-indent: 0pt)
      context {
        let n = ctr.get().first()
        text(weight: "bold", fill: accent)[#kind #n]
        if title != none { text(weight: "bold", fill: ieee_navy)[ (#title)] }
        text(weight: "bold", fill: ieee_navy)[.]
        if origin != none {
          h(4pt)
          _origin_badge(origin)
        }
        [ ]
      }
      emph(body)
    },
  )
}

#let theorem(title: none, origin: "Original Proof", body) = _stmt_box("Theorem", thm_counter, title, origin, body)
#let lemma(title: none, origin: "Original Proof", body) = _stmt_box("Lemma", lem_counter, title, origin, body)
#let proposition(title: none, origin: "Original Proof", body) = _stmt_box("Proposition", prop_counter, title, origin, body)
#let corollary(title: none, origin: "Original Proof", body) = _stmt_box("Corollary", cor_counter, title, origin, body)

#let proof(label: "Proof", body) = {
  block(
    width: 100%,
    above: 0.25em,
    below: 0.6em,
    breakable: true,
    {
      set par(first-line-indent: 0pt)
      text(style: "italic", weight: "bold", fill: ieee_slate)[#label: ]
      body
      h(1fr)
      text(fill: ieee_navy)[$square$]
    },
  )
}

#let remark(title: none, body) = {
  block(
    width: 100%,
    above: 0.45em,
    below: 0.5em,
    inset: (left: 6pt, top: 2pt, bottom: 2pt),
    stroke: (left: 1.2pt + ieee_border),
    breakable: true,
    {
      set par(first-line-indent: 0pt)
      text(weight: "bold", style: "italic", fill: ieee_slate)[Remark]
      if title != none { text(weight: "bold", style: "italic", fill: ieee_slate)[ (#title)] }
      text(weight: "bold", style: "italic", fill: ieee_slate)[: ]
      body
    },
  )
}

// =============================================================================
// IEEE Pseudocode Algorithm Box
// =============================================================================
#let algbox(title: "Algorithm", inputs: "", outputs: "", body) = {
  alg_counter.step()
  block(
    width: 100%,
    above: 0.6em,
    below: 0.65em,
    inset: (x: 5pt, y: 4.5pt),
    stroke: (top: 1.1pt + ieee_navy, bottom: 1.1pt + ieee_navy),
    fill: ieee_bg_alg,
    breakable: true,
    {
      set par(first-line-indent: 0pt, justify: false, leading: 0.46em, spacing: 0.34em)
      set text(size: 8.3pt)
      context {
        let n = alg_counter.get().first()
        text(weight: "bold", fill: ieee_navy)[Algorithm #n: #title]
      }
      v(2pt)
      line(length: 100%, stroke: 0.45pt + ieee_navy)
      v(2pt)
      if inputs != "" {
        [*Input:* #inputs \ ]
      }
      if outputs != "" {
        [*Output:* #outputs \ ]
      }
      v(2pt)
      body
    },
  )
}

// =============================================================================
// Full-width Figure / Wide Table Wrapper (Spans both IEEE columns cleanly)
// =============================================================================
#let wide_figure(caption: "", body) = {
  place(
    top + center,
    scope: "parent",
    float: true,
    block(
      width: 100%,
      above: 4pt,
      below: 8pt,
      {
        set align(center)
        body
        if caption != "" {
          v(3pt)
          set align(left)
          set par(first-line-indent: 0pt, justify: true)
          text(size: 8pt, fill: ieee_slate)[#caption]
        }
      },
    ),
  )
}
