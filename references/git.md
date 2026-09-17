# Git collaboration

Inspect state before changes. Commit or push only when requested. Exclude data, results, logs, checkpoints, caches, credentials, and private keys.

Prefer SSH. With old Git clients, verify identity with `ssh -T` and use a temporary `ssh-agent` when an explicitly named key is not loaded. Never print or commit private key contents.
