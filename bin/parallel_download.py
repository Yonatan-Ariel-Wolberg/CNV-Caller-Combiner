#!/usr/bin/env python3

import os
import sys
import urllib.request
import threading
import time

def download_segment(url, start_byte, end_byte, part_num, temp_dir):
    part_path = os.path.join(temp_dir, f"part_{part_num}")
    req = urllib.request.Request(url)
    req.add_header("Range", f"bytes={start_byte}-{end_byte}")
    
    attempts = 5
    for attempt in range(attempts):
        try:
            with urllib.request.urlopen(req, timeout=30) as response, open(part_path, "wb") as fh:
                fh.write(response.read())
            print(f"Segment {part_num} downloaded successfully.")
            return
        except Exception as e:
            print(f"Attempt {attempt+1} failed for segment {part_num}: {e}")
            time.sleep(2)
    raise Exception(f"Failed to download segment {part_num} after {attempts} attempts.")

def main():
    if len(sys.argv) < 3:
        print("Usage: parallel_download.py <URL> <Output_Path> [threads]")
        sys.exit(1)
        
    url = sys.argv[1]
    out_path = sys.argv[2]
    num_threads = int(sys.argv[3]) if len(sys.argv) > 3 else 16
    
    print(f"Retrieving headers from {url}...")
    req = urllib.request.Request(url, method="HEAD")
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            headers = resp.info()
            content_length = headers.get("Content-Length")
            accept_ranges = headers.get("Accept-Ranges")
    except Exception as e:
        print(f"Failed to get headers: {e}")
        sys.exit(1)
        
    if not content_length:
        print("Error: Server did not return Content-Length. Cannot download in parallel.")
        sys.exit(1)
        
    total_size = int(content_length)
    print(f"Total file size: {total_size / (1024*1024):.2f} MB")
    
    # If range requests are not accepted, fall back to single threaded
    if accept_ranges != "bytes":
        print("Warning: Server does not accept Accept-Ranges: bytes. Falling back to single-threaded download...")
        urllib.request.urlretrieve(url, out_path)
        print("Download complete.")
        sys.exit(0)
        
    temp_dir = out_path + ".parts"
    os.makedirs(temp_dir, exist_ok=True)
    
    segment_size = total_size // num_threads
    threads = []
    
    start_time = time.time()
    print(f"Launching {num_threads} threads for parallel download...")
    
    for i in range(num_threads):
        start_byte = i * segment_size
        # The last thread takes any remaining bytes
        end_byte = (i + 1) * segment_size - 1 if i < num_threads - 1 else total_size - 1
        
        t = threading.Thread(
            target=download_segment,
            args=(url, start_byte, end_byte, i, temp_dir)
        )
        threads.append(t)
        t.start()
        
    for t in threads:
        t.join()
        
    print("All segments downloaded. Stitching segments together...")
    with open(out_path, "wb") as final_fh:
        for i in range(num_threads):
            part_path = os.path.join(temp_dir, f"part_{i}")
            with open(part_path, "rb") as part_fh:
                final_fh.write(part_fh.read())
            os.remove(part_path)
            
    os.rmdir(temp_dir)
    
    print("Verifying archive integrity...")
    import tarfile
    try:
        with tarfile.open(out_path, "r:gz") as tar:
            # Attempt to read members list to verify Gzip block integrity
            tar.getnames()
        print("Archive integrity verified (not corrupted).")
    except Exception as e:
        print(f"Error: Downloaded archive is corrupted: {e}", file=sys.stderr)
        if os.path.exists(out_path):
            os.remove(out_path)
        sys.exit(2)
        
    duration = time.time() - start_time
    print(f"Download complete! Saved to {out_path} in {duration:.2f} seconds.")

if __name__ == "__main__":
    main()
