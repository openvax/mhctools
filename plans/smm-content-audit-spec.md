# SMM subset content audit (#465)

Keep the existing normalized GNU tar format and pin its uncompressed SHA-256.
The builder will verify that digest after writing the gzip archive, report both
checksums, and accept equivalent tar streams produced by different compressors.
Keep the official source checksum and the installer's exact hosted-artifact
checksum unchanged. No new manifest format or dependencies are needed.

Derive the content pin from the checksum-verified hosted artifact and confirm it
with fresh builds from the checksum-verified official archive on macOS and Linux.
Regression tests must accept different gzip encodings of the same tar stream,
reject content changes, and retain strict installer checksum validation. Remove
the builder's old placeholder checksum bypass and update the maintenance docs.

Release 4.4.3 after lint, full tests and green PR CI; deploy from clean master,
then verify the PyPI artifacts, release tag and deployment lock cleanup.
