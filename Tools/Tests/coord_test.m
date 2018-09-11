function pass = coord_test(trial_data)

    method = {'opensim','markers'};

    pass = false(4,1);

    for i = 1:2
        td_true = addSphereHand2TD(trial_data,method{i});
        td_test = addCoordPoint2TD(trial_data,struct('method',method{i},'point','hand','coord','sph'));

        truth = getSig(td_true,'sphere_hand_pos');
        test = getSig(td_test,'sph_hand_pos');

        pass(i) = all(all(truth == test));

        td_true = addCylHand2TD(trial_data,method{i});
        td_test = addCoordPoint2TD(trial_data,struct('method',method{i},'point','hand','coord','cyl'));

        truth = getSig(td_true,'cyl_hand_pos');
        test = getSig(td_test,'cyl_hand_pos');

        pass(2+i) = all(all(truth == test));
    end
