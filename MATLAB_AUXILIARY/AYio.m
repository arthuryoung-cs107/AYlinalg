classdef AYio
    properties (Constant)
        dat_dir_def = '../dat_dir/';
    end
    properties
        dat_dir = AYio.dat_dir_def;
    end
    methods
        function obj = AYio(dat_dir_)
            if (nargin == 1)
                dat_dir = dat_dir_;
            else
                dat_dir = AYio.dat_dir_def;
            end
        end
    end
    methods (Static)
        function write_matrix(mat_,name_,prefix_flag_)
            if (nargin==3)
                prefix_flag = prefix_flag_;
            else
                prefix_flag = true;
            end
            if (prefix_flag)
                name_full = [AYio.dat_dir_def name_ '.aydat'];
            else
                name_full = [name_ '.aydat'];
            end

            mat_size = size(mat_);
            header = [2 2 prod(mat_size)];

            file_id = fopen(name_full,'w+');
            fwrite(file_id,header,'int');
            fwrite(file_id,mat_size,'int');
            fwrite(file_id,reshape(mat_',[],1),'double');
            fclose(file_id);
        end
        function [mat_out out1 out2] = read_matrix(name_,prefix_flag_)
            if (nargin==2)
                prefix_flag = prefix_flag_;
            else
                prefix_flag = true;
            end
            if (prefix_flag)
                name_full = [AYio.dat_dir_def name_ '.aydat'];
            else
                name_full = [name_ '.aydat'];
            end

            file_id = fopen(name_full);
            header = fread(file_id,5,'int=>int');
            dims = header(4:end);
            mat_raw = fread(file_id,prod(dims),'double=>double');
            fclose(file_id);

            mat_out = (reshape(mat_raw,dims(2),dims(1)))';

            if (nargout>1)
                if (nargout==2)
                    out1 = dims;
                else
                    [out1 out2] = deal(dims(1),dims(2));
                end
            end
        end
    end
end
