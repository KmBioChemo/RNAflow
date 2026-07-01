# Linked interactive volcano

A crosstalk-linked volcano plot: brushing points in the plotly volcano
highlights the matching rows in a DT table (and vice versa), with the
selection also surfaced server-side. The data-prep is a pure, testable
function; the figure builder wraps a
[`crosstalk::SharedData`](https://rdrr.io/pkg/crosstalk/man/SharedData.html).
