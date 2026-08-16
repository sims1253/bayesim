# Internal lifecycle contracts

The public simulation API intentionally keeps the historical `skipped`
status for compatibility. A skipped task is not, however, a completed
scientific outcome: `stop_reason` records why execution stopped and
resume turns policy-stopped tasks back into pending work.
