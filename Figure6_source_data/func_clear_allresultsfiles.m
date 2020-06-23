function func_clear_oldfiles()

    par = {'mouse1_session1';
           'mouse1_session2';
           'mouse1_session3';
           'mouse1_session4';

           'mouse2_session1';
           'mouse2_session2';
           'mouse2_session3';
           'mouse2_session4';
           };

    stepflist = [5:15];

    for ip = [1,2,3,4,5,6,7,8]
        session = par{ip,1};
        for istepf = 1:length(stepflist)
            disp(['running session: ',session,' ...']);
            delete(sprintf('%s\\*.fig',session));
%             delete(sprintf('%s\\correct_DecodingPrediction*',session));
%             delete(sprintf('%s\\incorrect_DecodingPrediction*',session));
%             delete(sprintf('%s\\correct_SVMClassification*',session));
%             delete(sprintf('%s\\incorrect_SVMClassification*',session));
        end
    end
end

