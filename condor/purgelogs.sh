#!/bin/bash

# Clean log/errors/output of junk
# Some functionality in the future - check contents of a folder, then if empty, return "folder is empty"
# vs. if full, would you like to delete files?

# clear log contents
echo "Clearing log folder..."
rm log/*.log

# clear errors contents
echo "Clearing errors folder..."
rm errors/*.err

# clear output contents
echo "Clearing output folder..."
rm output/*.out

# clear Results contents
echo "Clearing Results folder..."
rm Results/*.csv

echo ""; echo "...done."; echo ""; echo "Relevant directories have been emptied:"

./checkstuff.sh
