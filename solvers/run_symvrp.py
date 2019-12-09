problem_file = argv[1]
time_limit = 10

# plan A: get from filename
if 

# plan B: get from VEHICLES field of the file


# plan C: get using the capacity

m = regexp_k_from_filename.search( path.basename(vrp_file) )
        k = None
        if m:
            k = int(m.group(1))
        elif K:
            k = K
        else:
            k = ceil(sum(demands)*1.05/C) # estimate k by giving 5 % slack capacity
	    

with open(problem_file, 'r'

cmd_w_args = ["vrp.exe",
          '-F', problem_file,
          '-p', '1',
          '-N', str(k), 
          '-t', str(time_limit),
          '-v', str(verbosity)
    ]
    
if upper_bound:
	ub_at_end = upper_bound
	cmd_w_args.append('-u')
	cmd_w_args.append(str(upper_bound+0.1))

p = Popen(cmd_w_args, stdout=PIPE, bufsize=1)

