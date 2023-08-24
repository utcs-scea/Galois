#####   Inputs   #####
# hosts = [8, 16, 32, 64]
# scales = ["31", "31", "31", "31"]
hosts = [4]
scales = ["5"]



#####   Parameters   #####
data_size = 8
id_size = 8
addr_size = 8



#####   Processing   #####
for host in hosts:
	host_index = hosts.index(host)
	scale = scales[host_index]

	local_read = []
	master_read = []
	master_write = []
	mirror_read = []
	mirror_write = []
	dirty_mirror_to_remote = []
	
	for host_id in range(host):
		local_read.append([])
		master_read.append([])
		master_write.append([])
		mirror_read.append([])
		mirror_write.append([])
		dirty_mirror_to_remote.append([])
		
		round_curr = -1
		round_total = 0
		
		filename = "bfs_" + scale + "_" + str(host) + "procs_id" + str(host_id)
		f = open(filename, "r")
		lines = f.readlines()
		
		for row in lines:
			if row.find('#####   Round') != -1:
				row_split = row.strip().split()
				round_curr += 1
				round_total += 1

				local_read[host_id].append([])
				master_read[host_id].append([])
				master_write[host_id].append([])
				mirror_read[host_id].append([])
				mirror_write[host_id].append([])
				dirty_mirror_to_remote[host_id].append([-1 for i in range(host)])
			elif row.find('local read (stream)') != -1:
				row_split = row.strip().split()
				num = row_split[5]
				local_read[host_id][round_curr] = int(num)
			elif row.find('master reads') != -1:
				row_split = row.strip().split()
				num = row_split[4]
				master_read[host_id][round_curr] = int(num)
			elif row.find('master writes') != -1:
				row_split = row.strip().split()
				num = row_split[4]
				master_write[host_id][round_curr] = int(num)
			elif row.find('mirror reads') != -1:
				row_split = row.strip().split()
				num = row_split[4]
				mirror_read[host_id][round_curr] = int(num)
			elif row.find('mirror writes') != -1:
				row_split = row.strip().split()
				num = row_split[4]
				mirror_write[host_id][round_curr] = int(num)
			elif row.find('remote communication for') != -1:
				row_split = row.strip().split()
				to_host = int(row_split[6][:-1])
				num = row_split[7]
				dirty_mirror_to_remote[host_id][round_curr][to_host] = int(num)
		
		f.close()

	filename = "GAL_BFS_" + str(host) + "_mirror100_" + scale + ".stats"
	f = open(filename, "w")

	for round_num in range(round_total):
		f.write("BSP " + str(2*round_num) + "\n")

		for src in range(host):
			for dst in range(host):
				if (src == dst): # only local in compute phase`
					stream_read = local_read[src][round_num]
					stream_read_bytes = stream_read * addr_size
					
					random_read = master_read[src][round_num] + mirror_read[src][round_num]
					random_read_bytes = random_read * data_size
					
					random_write = master_write[src][round_num] + mirror_write[src][round_num]
					random_write_bytes = random_write * data_size
					
					if stream_read != 0 or random_read != 0 or random_write != 0:
						f.write("# " + str(src) + " " + str(dst) + "\n")
					
					if stream_read != 0:
						f.write("STR RD " + str(stream_read) + " " + str(stream_read_bytes) + "\n")

					if random_read != 0:
						f.write("RND RD " + str(random_read) + " " + str(random_read_bytes) + "\n")
					
					if random_write != 0:
						f.write("RND WR " + str(random_write) + " " + str(random_write_bytes) + "\n")
		
		f.write("\n")

		f.write("BSP " + str(2*round_num+1) + "\n")

		for src in range(host):
			for dst in range(host):
				if (src != dst): # only remote in communication phase
					random_write = dirty_mirror_to_remote[src][round_num][dst]
					random_write_bytes = random_write * data_size
					
					if random_write != 0:
						f.write("# " + str(src) + " " + str(dst) + "\n")
						f.write("RND RMW " + str(random_write) + " " + str(random_write_bytes) + "\n")

		f.write("\n")
	
	f.close()
