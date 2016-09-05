
-- Create software HUMAnN2, the humann2 command and it's parameters
-- magic # 1 because this is an "artifact transformation" plugin

-- taken from:
-- http://stackoverflow.com/q/21386772

WITH row_software AS (
    INSERT INTO qiita.software (name, version, description, environment_script, start_script, software_type_id) 
    VALUES ('HUMAnN2', '0.9.1', 'HUMAnN2 is the next generation of HUMAnN (HMP Unified Metabolic Analysis Network)', 'workon qp_shotgun', 'start_qp_shotgun', 1) 
    RETURNING software_id AS sid
), row_software_command AS (
    INSERT INTO qiita.software_command (software_id, name, description) 
    VALUES ((SELECT sid FROM row_software), 'HUMAnN2 pipeline', 'HUMAnN2 pipeline') 
    RETURNING command_id AS cid
) 
INSERT INTO qiita.command_parameter (command_id, parameter_name, parameter_type, required, default_value) VALUES
((SELECT cid FROM row_software_command), 'input_data', 'artifact', True, NULL),
((SELECT cid FROM row_software_command), '--protein-database', 'string', True, 'uniref'),
((SELECT cid FROM row_software_command), '--nucleotide-database', 'string', True, 'chocophlan');

-- Inserting workflow
INSERT INTO qiita.default_workflow (software_id, name)
    VALUES ((SELECT software_id FROM qiita.software WHERE name='HUMAnN2' and version='0.9.1' ORDER BY software_id DESC LIMIT 1), 'HUMAnN2 workflow');

-- Inserting workflow default parameter set
INSERT INTO qiita.default_parameter_set (command_id, parameter_set_name, parameter_set)
    VALUES ((SELECT command_id FROM qiita.software_command WHERE name='HUMAnN2 pipeline' ORDER BY command_id DESC LIMIT 1), 'HUMAnN2 workflow - defaults', ('{"--protein-database": "uniref", "--nucleotide-database": "chocophlan"}')::json);

-- On later day we will need to add nodes, edges and their connections to the workflows