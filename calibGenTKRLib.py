# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['calibGenTKR'])
    env.Tool('addLibrary', library = env['rootGuiLibs'])
def exists(env):
    return 1;
