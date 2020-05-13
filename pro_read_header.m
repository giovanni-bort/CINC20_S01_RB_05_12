% function pro_read_header



	[H_recording,Total_time,H_num_leads,H_Fs,H_gain,H_age,H_sex]=extract_data_from_header(Header);

    fprintf('recording,time');fprintf(' %s %10.2f ',H_recording,Total_time);
    fprintf('leads:%6.0f Fs:%8.1fage:sex%6.0f%6.0f\n',H_num_leads,H_Fs,H_age,H_sex);
    fprintf('Gain:');fprintf('%8.1f',H_gain);fprintf('\n');
    