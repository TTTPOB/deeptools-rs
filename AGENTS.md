## ref impl test plan
- use pixi to create python test environment
- include bioconda channel, deeptools 3.5.6 package
- find relavant test data, generate an reference output using `pixi run computeMatrix ...`

## workspace
for openai codex: you are configured to run code in sandbox so that network access may not be reliable. when you try to run `cargo update` or `cargo metadata`, you may see errors. for `cargo update`, you don't need to run it. `rust-analyzer` will handle by automatically lock. `cargo metadata`, you may want to use it to get doc or codes, try to use rg to find them under local directory, as `rust-analyzer` may have already indexed them. 