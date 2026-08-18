# Internal run-store seam

The execution engine currently uses the established checkpoint
functions, but all durable state can be addressed through this small
adapter. Keeping this seam separate makes lifecycle tests independent of
the filesystem and gives future append-only stores a single contract to
implement.
