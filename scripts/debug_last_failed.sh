#!/usr/bin/env bash

# Helper script to automatically find the most recent failed Nextflow task directory,
# extract .command.run and .command.sh files, and prep them for local execution in WSL.

set -euo pipefail

WORK_DIR="work"

if [ ! -d "$WORK_DIR" ]; then
    echo "Error: Nextflow 'work' directory not found in the current folder."
    exit 1
fi

echo "Searching for the most recent failed task..."

# Find all .exitcode files, filter for non-zero values, and sort by modification time (newest first)
failed_exitcode_file=$(find "$WORK_DIR" -name ".exitcode" -type f | while read -r f; do
    code=$(cat "$f" 2>/dev/null || echo "0")
    if [ "$code" -ne 0 ]; then
        echo "$(stat -c %Y "$f") $f"
    fi
done | sort -nr | head -n1 | cut -d' ' -f2-)

if [ -z "$failed_exitcode_file" ]; then
    echo "No failed tasks found in '$WORK_DIR/' (all exit codes are 0 or no tasks have run)."
    exit 0
fi

task_dir=$(dirname "$failed_exitcode_file")
echo "Found failed task directory: $task_dir"
echo "Exit code: $(cat "$failed_exitcode_file")"

# Create a debug directory
debug_dest="debug_task"
mkdir -p "$debug_dest"
echo "Extracting files to './$debug_dest/'..."

# Copy command files
for file in .command.sh .command.run .command.err .command.out .command.log .exitcode; do
    if [ -f "$task_dir/$file" ]; then
        cp "$task_dir/$file" "$debug_dest/$file"
        # Ensure scripts are executable
        if [[ "$file" == *".sh" || "$file" == *".run" ]]; then
            chmod +x "$debug_dest/$file" 2>/dev/null || true
        fi
        # Convert CRLF to LF for WSL compatibility
        sed -i 's/\r$//' "$debug_dest/$file" 2>/dev/null || true
    fi
done

echo "--------------------------------------------------"
echo "Preparation complete! Files copied to './$debug_dest/'"
echo "To run the task locally in WSL, you can execute:"
echo "  cd $debug_dest && ./.command.run"
echo "Or run the raw commands directly via:"
echo "  cd $debug_dest && ./.command.sh"
echo "--------------------------------------------------"
