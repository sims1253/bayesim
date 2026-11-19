# Internal run-store seam

The run store keeps the latest checkpoint state between writes and uses
the checkpoint functions to read and persist outcomes. Runs without a
result path keep outcomes in the execution loop and do not use a store.
