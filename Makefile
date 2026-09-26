# Convenience wrapper around MATLAB's buildtool (see buildfile.m).
# Windows users can run the same tasks from MATLAB: buildtool <task>.
MATLAB ?= matlab

.PHONY: test docs examples manual release clean

test docs examples manual release clean:
	$(MATLAB) -batch "buildtool $@"
