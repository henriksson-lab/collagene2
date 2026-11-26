# Native saved-catalog backup controls

`CatalogBackupControls` is available only in the native application. Root injects a typed native callback; the view never handles backup archive bytes, browser downloads, HTTP, operator paths or arbitrary file names. It has its own busy state and invalidates pending review/result state when its local lab context or lock changes, or the component closes.

The scope is the saved local catalog: plasmids, folders, saved revisions and durable repository state. Users must save work in the local library first. These controls do not back up the current autosaved working-editor files. Creating a backup does not switch or alter the current repository.

Known backups are shown by their managed receipt basename, UTC creation date and archive size, with storage counts and archive/database/schema hashes under receipt details. Pagination retains the native catalog revision; refresh resets the snapshot. A known unfinished creation intent remains visible without being described as a completed backup. Users can request native verification of its exact known identity; only completed known receipts supply an expected archive hash.

A list or creation receipt does not enable restore. A matching native Verify result creates a local restore review pinned to the known backup ID, basename, hash and complete local lab context. The user must separately confirm restoration to a new managed copy. The backend rechecks the archive and writes that distinct copy; the UI requires its real matching recovery receipt before reporting success. It never overwrites, merges or switches the current catalog or replaces the working editor. Storage-integrity verification is not molecular-design or experimental verification.

The pure `catalog_backup_model::tests` fixtures cover completed-entry hash binding, unfinished-intent verification without direct restoration, local-scope changes and substituted recovery receipts. The view itself is native/WASM frontend code; no browser automation is required by these fixtures.
