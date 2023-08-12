all: clippy test

clippy:
	cargo clippy

test:
	cargo test
