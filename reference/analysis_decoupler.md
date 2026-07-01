# Transcription-factor and pathway activity inference

Infers per-sample-independent transcription-factor (CollecTRI regulons)
and pathway (PROGENy) *activity* scores from a differential-expression
contrast, using decoupleR. Instead of asking "which genes changed",
activity inference asks "which upstream regulators / pathways best
explain the change", by scoring a prior-knowledge network against the
ranked DE statistic (a univariate linear model for TFs, a multivariate
one for pathways).

## Details

The scoring is a pure function of the DE table and the network; only the
network *fetch*
([`get_tf_network()`](https://KmBioChemo.github.io/RNAflow/reference/get_tf_network.md)
/
[`get_pathway_network()`](https://KmBioChemo.github.io/RNAflow/reference/get_pathway_network.md))
reaches OmniPath over the internet.
