# Convenience wrapper around MATLAB's buildtool (see buildfile.m).
# Windows users can run the same tasks from MATLAB: buildtool <task>.
MATLAB ?= matlab

.PHONY: test docs examples manual release package clean

test docs examples manual release package clean:
	$(MATLAB) -batch "buildtool $@"
