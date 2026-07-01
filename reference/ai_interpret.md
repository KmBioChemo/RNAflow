# AI-assisted biological interpretation

Turns a differential-expression contrast (plus optional functional
enrichment) into a compact prompt and asks Anthropic's Claude API to
write a biological narrative. The prompt-building functions are pure and
testable without any network access;
[`call_claude()`](https://KmBioChemo.github.io/RNAflow/reference/call_claude.md)
is the only function that touches the API (guarded by httr2).

## Details

Only a compact summary leaves the machine – top gene *names* with their
fold-changes and FDRs, and the top enrichment terms. Count matrices and
sample metadata are never sent.
