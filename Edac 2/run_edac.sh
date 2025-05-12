#!/bin/bash
# Run edac with input file and kill it after "That's all, folks!" appears
./edac < edac.in > edac_output.tmp 2>&1 &
pid=$!

# Wait for "That's all, folks!" to appear in the output
while true; do
  if grep -q "That's all, folks!" edac_output.tmp; then
    # Message found, kill the process
    kill $pid 2>/dev/null
    wait $pid 2>/dev/null
    break
  fi
  
  # Check if process is still running
  if ! ps -p $pid > /dev/null; then
    break
  fi
  
  # Sleep for a short time before checking again
  sleep 0.5
done

# Display the output
cat edac_output.tmp

# Clean up
rm edac_output.tmp 