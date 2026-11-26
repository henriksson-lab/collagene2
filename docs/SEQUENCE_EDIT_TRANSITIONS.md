# Sequence editor navigation after edits

Keyboard replacement, paste and cut commit through the active document’s ordinary history controller. After a successful commit, the sequence view places its caret using the change in normalized sequence length. Raw clipboard bytes do not determine a coordinate. The cursor records the new source event and revision after that synchronous commit.

App’s later source invalidation effect distinguishes these local edits from external changes. Its private navigation tracker binds the workspace epoch, native document ID, exact head event, immutable revision and sequence length before and after the commit. Only contiguous local transitions rooted at the previously observed head preserve the resulting caret. Multiple synchronous edits can form one exact chain before the next effect tick. This is local presentation state; it creates no scientific or physical authority and is never serialized into an archive.

A changed active document, external edit, equal-content event replacement or restored workspace clears selection. Successful undo/redo also clears the pending transition and selection. A rejected history operation leaves the current caret alone. Writes to other documents leave active selection intact while still invalidating multi-source plan reports. The tracker starts with the actual active head, so the first effect tick does not clear a newly placed pointer selection. A pointer selection made after a local edit takes precedence over the earlier caret.

The source effect reads caret state without tracking it, and uses the immutable head’s cached revision instead of hashing every open molecule. Sequence view identity includes the exact head event, so equal-content events cannot retain an old keyboard anchor.

`sequence_edit_transition.rs` contains six native fixture groups using actual `DocumentHistory` replacement, commit, fork, undo and redo operations. They cover repeated characters, coalesced edits, whitespace-normalized paste, cut, no-op replacement, unrelated document writes, equal-content events/documents, restore epochs, external interleaving and rejected geometry. Browser checks are not part of this change.
