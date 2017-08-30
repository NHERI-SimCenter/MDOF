classdef ElasticPPSpring < handle

    % an ElasticPP spring
    % written: fmk 10/2016
    properties
        Ke;   % Elastic tangent
        Fy;   % Yield Force
        eyPos; % displacement Marking Pos Yield for last committed step
        eyNeg; % displacement MArking Neg Yield for last committed step
        currentDeformation; % current Trial Deformation
    end
    
    methods
        % constructor
        function obj =  ElasticPPSpring(ke,fy)
            obj.Ke = ke;
            obj.Fy = fy;
            obj.eyPos = fy/ke;
            obj.eyNeg = -obj.eyPos;
            obj.currentDeformation = 0;
        end
            
        % method to set current trial deformation & obtain K and fs
        function [K, fs] = setTrialDisplacement(obj, deformation)
            obj.currentDeformation = deformation;
            if deformation > obj.eyPos
                K = 0;
                fs = obj.Fy;
            elseif deformation < obj.eyNeg
                K = 0;
                fs = -obj.Fy;
            else
                K = obj.Ke;
                fs = -obj.Fy + (deformation-obj.eyNeg)*obj.Ke; 
            end
        end
        
        %method to commitState, i.e. accepted solution on path
        function commitState(obj)
            if obj.currentDeformation > obj.eyPos
                obj.eyPos = obj.currentDeformation;
                obj.eyNeg = obj.eyPos - 2*obj.Fy/obj.Ke;
            end
            if obj.currentDeformation < obj.eyNeg
                obj.eyNeg = obj.currentDeformation;
                obj.eyPos = obj.eyNeg + 2*obj.Fy/obj.Ke;
            end 
        end
    end
    
end

