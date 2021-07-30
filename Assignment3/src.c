#include<stdio.h>
#include<mpi.h>
#include<stdlib.h>

#define FLOAT_MAX 3.40282e+038

FILE* file;
FILE* out;
FILE* dump;

// returns number of columns in csv file
long long get_cols_num(){
	int c;
	long long cols = 1;
	c = fgetc(file);
	while (c != '\n'){
		if (c==',') cols++;
		c = fgetc(file);
	}
	return cols;
}

// returns number of rows (apart from header) in s=csv file
long long get_rows_num(){
	char c;
	long long rows = -1;
	for (char c=fgetc(file); c!=EOF; c=fgetc(file)){
		if (c == '\n') rows++;
	}
	fseek(file, 0, SEEK_SET);
	return rows;
}

int main(int argc, char* argv[]){

	MPI_Init(NULL,NULL);

	long long block_size; // stores the block size to be communicated and computed at a time
	long long rows, cols; // rows and cols in csv file

	int size, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	float* data;
	if (rank == 0){
		file = fopen(argv[1], "r");
		out = fopen("output.txt", "w"); // this file shall have the output as required
		dump = fopen("dump", "a"); // this file shall have raw time readings for plotting purpose

		rows = get_rows_num();
		cols = get_cols_num();
		cols = cols - 2; // to remove first two columns (latitude and longitude)
		block_size = 1000; // optimal block size (based on experiments)

		data = (float*)malloc(sizeof(float*)*rows*cols);

		// scan the entire data in 1D array "data"
		for (long long i=0; i<rows; i++){
			float dump;
			fscanf(file, "%f,%f,", &dump, &dump);
			for (long long j=0; j<cols-1; j++){
				fscanf(file, "%f,", &data[i*cols + j]);
			}
			fscanf(file, "%f\n", &data[i*cols + cols-1]);
		}
	}

	// if there is a single process, no parallelization possible
	if (size == 1){
		double s_time = MPI_Wtime();

		float min_val[cols]; // will store yearwise minimum values
		for (long long i=0; i<cols; i++) min_val[i] = FLOAT_MAX;

		float global_min = FLOAT_MAX; // will store global minimum
		for (long long i=0; i<rows; i++){
			for (long long j=0; j<cols; j++){
				min_val[j] = min_val[j] > data[i*cols + j] ? data[i*cols + j] : min_val[j];
			}
		}
		for (long long i=0; i<cols; i++){
			global_min = global_min > min_val[i] ? min_val[i] : global_min;
		}

		double e_time = MPI_Wtime();

		// output required information to output.txt
		fprintf(out, "%f", min_val[0]);
		for (long long i=1; i<cols; i++){
			fprintf (out, ",%f", min_val[i]);
		}
		fprintf (out, "\n%f\n%f", global_min, e_time-s_time);
		fprintf (dump, "%f\n", e_time-s_time);
	}

	// if multiple processes spawned, parallelization possible using 1D decomposition along rows
	if (size > 1){
		// if barrier is not put then other processes will keep waiting for root process which mght be still reading data
		MPI_Barrier(MPI_COMM_WORLD);
		double s_time = MPI_Wtime();

		long long meta_data[3];
		if (!rank){
			meta_data[0] = block_size;
			meta_data[1] = rows;
			meta_data[2] = cols;
		}
		MPI_Bcast(meta_data, 3, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
		if (rank){
			block_size = meta_data[0];
			rows = meta_data[1];
			cols = meta_data[2];
		}
		
		long long blocks = rows/block_size; // total blocks of data with block_size being the number of rows in each block
		long long rem_data_offset = blocks*block_size; // offset for block of data that has rows < block_size
		long long blocks_per_process = blocks/size; // equal distribution of all the blocks

		// one buffer shall be receiving data while the other buffer is used for computation
		float* recv_data[2];
		recv_data[0] = (float*)malloc(sizeof(float*)*block_size*cols);
		recv_data[1] = (float*)malloc(sizeof(float*)*block_size*cols);
		
		if (!rank){
			float min_val[cols]; // will store local min vals (calculated by root on its share of data)
			for (long long i=0; i<cols; i++) min_val[i] = FLOAT_MAX;

			/* 
			data distribution strategy ->
			scatter 1st block of data
			loop 'blocks_per_process' number of times :
				issue an iscatter to distribute (n+1)th block
				perform min val computation on (n)th block, which is already received
				issue mpi_wait to ensure receival of (n+1)th block
			*/
			MPI_Scatter(data, block_size*cols, MPI_FLOAT, recv_data[0], block_size*cols, MPI_FLOAT, 0, MPI_COMM_WORLD);
			for (long long i=1; i<blocks_per_process; i++){
				MPI_Request req;
				MPI_Iscatter(data + i*block_size*cols, block_size*cols, MPI_FLOAT, recv_data[i%2], block_size*cols, MPI_FLOAT, 0, MPI_COMM_WORLD, &req);
				long long offset = (i-1)*block_size*size;
				for (long long j=offset; j<offset + block_size; j++){
					for (long long k=0; k<cols; k++){
						min_val[k] = min_val[k] > data[j*cols + k] ? data[j*cols + k] : min_val[k];
					}
				}
				MPI_Wait(&req, MPI_STATUS_IGNORE);
			}
			long long offset = (blocks_per_process-1)*block_size*size;
			for (long long j=offset; j<offset + block_size; j++){
				for (long long k=0; k<cols; k++){
					min_val[k] = min_val[k] > data[j*cols + k] ? data[j*cols + k] : min_val[k];
				}
			}
			// remaining data calculation performed on root since distributing even this will lead to communication overhead
			for (long long j=rem_data_offset; j<rows; j++){
				for (long long k=0; k<cols; k++){
					min_val[k] = min_val[k] > data[j*cols + k] ? data[j*cols + k] : min_val[k];
				}
			}
			// mpi reduce for final mins calculation (receives local mins calculated by every process)
			float final_mins[cols];
			MPI_Reduce(min_val, final_mins, cols, MPI_FLOAT, MPI_MIN, 0, MPI_COMM_WORLD);

			// global min calculation performed only on root as num of cols <= 100
			float global_min = FLOAT_MAX;
			for (long long i=0; i<cols; i++){
				if(global_min > final_mins[i]) global_min = final_mins[i];
			}

			double time = MPI_Wtime() - s_time;

			// mpi reduce to get max times
			double max_time;
			MPI_Reduce(&time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

			// output required information to output.txt
			fprintf(out, "%f", final_mins[0]);
			for (long long i=1; i<cols; i++){
				fprintf (out, ",%f", final_mins[i]);
			}
			fprintf (out, "\n%f\n%f", global_min, max_time);
			fprintf (dump, "%f\n", max_time);
		}
		if (rank){
			MPI_Scatter(data, block_size*cols, MPI_FLOAT, recv_data[0], block_size*cols, MPI_FLOAT, 0, MPI_COMM_WORLD);
			float min_val[cols]; // will store local min vals (calculated by this process on its share of data)
			for (long long i=0; i<cols; i++) min_val[i] = FLOAT_MAX;
			
			for(long long i=1; i<blocks_per_process; i++){
				MPI_Request req;
				MPI_Iscatter(data + i*block_size*cols, block_size*cols, MPI_FLOAT, recv_data[i%2], block_size*cols, MPI_FLOAT, 0, MPI_COMM_WORLD, &req);
				int data_label = (i-1)%2;
				for (long long j=0; j<block_size; j++){
					for (long long k=0; k<cols; k++){
						min_val[k] = min_val[k] > recv_data[data_label][j*cols + k] ? recv_data[data_label][j*cols + k] : min_val[k];
					}
				}
				MPI_Wait(&req, MPI_STATUS_IGNORE);
			}
			int data_label = (blocks_per_process-1)%2;
			for (long long j=0; j<block_size; j++){
				for (long long k=0; k<cols; k++){
					min_val[k] = min_val[k] > recv_data[data_label][j*cols + k] ? recv_data[data_label][j*cols + k] : min_val[k];
				}
			}
			// mpi reduce for final mins calculation
			double final_mins[cols];
			MPI_Reduce(min_val, final_mins, cols, MPI_FLOAT, MPI_MIN, 0, MPI_COMM_WORLD);

			double time = MPI_Wtime() - s_time;
			// mpi reduce to get max times
			double max_time;
			MPI_Reduce(&time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		}
	}
	MPI_Finalize();
	return 0;
}