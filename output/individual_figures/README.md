Every file in this folder, `individual_figures`, is handwritten.
In most cases, they're table headers and notes, meant to be used as a wrapper around the table bodies in `tex_fragments`.

For example, the file `summary_stats.tex` is meant to be used like the following.
Note that there's no `\caption` or `\label`; these are already included in `summary_stats.tex`.
In turn, `summary_stats.tex` includes fragments for each panel.
Note that `summary_stats.tex` should use the command `\input` and needs to specify paths *relative to itself*, not relative to `all_figures.tex`.
(I'm unclear on why the `\subimport` command doesn't work here.)

Dummy contents of `main_document.tex`:
```latex
...
\begin{table}[!hbt]
\import{individual_figures/}{summary_stats.tex}
\end{table}
```

Dummy contents of `summary_stats.tex`:
```latex
...
\begin{table}[!hbt]
\begin{tabular}{lll}
\input{../tex_fragments/summary_stats_panel_A.tex}
\input{../tex_fragments/summary_stats_panel_B.tex}
\end{tabular}
\end{table}
```
