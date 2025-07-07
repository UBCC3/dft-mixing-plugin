## Database

- Resolve coefficeints as part of preprocessing, so
PSI4 does not really have to do the heavy lifting resolving the coefficeints (and also more efficient retrieval from the database.)


## Computation Side
- Implement LibXC functional decomposition, to save more computation time. This will need a heavy overhaul on the Python side especially for base functionals, need to make it much more cleaner, as base functional resolution is still rather weird...

(I need the tweaks)

- Test LCOM implementation
