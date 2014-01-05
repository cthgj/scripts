


<!DOCTYPE html>
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en"
    class="
    
    "

>
<head>
    <meta http-equiv="X-UA-Compatible" content="IE=9" >
    <script type="text/javascript">window.ST=new Date().getTime();</script>
    <meta content="text/html; charset=UTF-8" http-equiv="content-type" />
    <meta name="description" content="Dropbox is a free service that lets you bring your photos, docs, and videos anywhere and share them easily. Never email yourself a file again!"/>
    <meta name="keywords" content="online storage, free storage, file sharing, share files, awesome, cloud storage, online backup, cross platform, sync, sharing, mac, windows, os x, linux, backup, collaboration, file versioning, file revisions, remote access, undelete"/>
    <meta name="google-site-verification" content="TnuSyOnBMNmtugbpL1ZvW2PbSF9LKvoTzrvOGS9h-b0" />
    <meta name="google-site-verification" content="EZKIczQcM1-DVUMz8heu1dIhNtxNbLqbaA9-HbOnCQ4" />
    <meta name="norton-safeweb-site-verification" content="tz8iotmk-pkhui406y41y5bfmfxdwmaa4a-yc0hm6r0fga7s6j0j27qmgqkmc7oovihzghbzhbdjk-uiyrz438nxsjdbj3fggwgl8oq2nf4ko8gi7j4z7t78kegbidl4" />
    <meta name="robots" content="noindex">
    <title>
                stats_for_precision.py -
        Dropbox
    </title>

    <link rel="shortcut icon" href="/static/26613/images/favicon.ico"/>

    <link href="/static/26613/css/main.css" rel="stylesheet" type="text/css"/>

    <link rel="apple-touch-icon" href="/static/26613/images/dropbox_webclip.png"/>
    <link rel="P3Pv1" href="/w3c/p3p.xml"/>



<script type="text/javascript">
var Constants = {
  BLOCK_CLUSTER: 'dl\x2dweb\x2edropbox\x2ecom',
  PUBSERVER: 'dl\x2edropbox\x2ecom',
  WEBSERVER: 'www\x2edropbox\x2ecom',
  LIVE_TRANSCODE_SERVER: 'showbox\x2dtr\x2edropbox\x2ecom',
  DISABLE_VIDEO_ICONS: false,
  DISABLE_VIDEOS_IN_LIGHTBOX: false,
  block: 'dl\x2dweb\x2edropbox\x2ecom',
  protocol: 'https',
  uid: '14062712',
  email: 'selenop\x40gmail\x2ecom',
  sess_id: '96404536493528980311302326560028675281',
  root_ns: 21921772,
  SVN_REV: '26613',
  TOKEN: 'daf38c04e4',
  IS_PROD: 1,
  WIT_ENABLED: 0,
  upload_debug: false,
  tcn: 'touch',
  date_format: 'M\x2fd\x2fyyyy',
  time_format: 'h\x3amm\x20a',
  datetime_format: 'M\x2fd\x2fyyyy\x20h\x3amm\x20a',
  ADMIN: 0,
  referrer: 'https\x3a\x2f\x2fwww\x2edropbox\x2ecom\x2fsh\x2fq8v738nqie29asu\x2fx9b6fqLDX4', 
  TWO_ITEM_LIST: '\x25\x28x\x29s\x20and\x20\x25\x28y\x29s',
  THREE_ITEM_LIST: '\x25\x28x\x29s\x2c\x20\x25\x28y\x29s\x2c\x20and\x20\x25\x28z\x29s',
  LOCALES: [["en", "English"], ["es", "Espa\u00f1ol"], ["fr", "Fran\u00e7ais"], ["de", "Deutsch"], ["ja", "\u65e5\u672c\u8a9e"]],
  USER_LOCALE: 'en',
  EMAIL_VERIFIED: 1
  ,
  auto_timezone_offset: 2.0
    ,
    can_shmodel: true
    ,
    can_fb_invite: true
};
</script>
    <script type="text/javascript">
        function global_report_exception (e, f, l, tb, force) {
            if (!window.reported_exception || force) {
                var stack_str = "";
                try {
                    if (!tb) {
                        var stack = get_stack_rep();
                        stack.pop(); 
                        stack.pop(); 
                        stack_str = stack.join("\n");
                    }
                } catch (e) { }
                var dbr = new Ajax.Request("/jse", {parameters: {'e': e, 'f': f || window.location.href, 'l': l, 'loc': window.location.href, 'ref': Constants.referrer, 'tb': tb || stack_str, 'trace': Trace && Trace.get() }});
                window.reported_exception = true;
           }
        }
        window.LoadedJsSuccessfully = false; 

        window.onerror = function (e, f, l) {
            global_report_exception(e, f, l);
        };
    </script>
<script type="text/javascript" src="/static/26613/javascript/dropbox-mini.js"></script>

    <script type="text/javascript">
      if (self != top) {
          top.location.replace(self.location.href);
          setTimeout(function() {document.body.innerHTML = "<img src='https://www.dropbox.com/static/images/logo.png' onClick='top.location.href=window.location.href'>";},1);
      }
    </script>

    <style type="text/css">
      .hny-yvsukf { display: none; }
    </style>
    <!--[if lte IE 8]> 
    <style>
    .sick-input input[type=password] {
      padding-top: 7px !important;
      padding-bottom: 3px !important;
    }
    .sick-input.small input {
      padding-top: 4px !important;
      padding-bottom: 2px !important;
    }
    .sick-input.small input[type=password] {
      padding-top: 6px !important;
      padding-bottom: 0px !important;
    }
    </style>
    <![endif]-->
    <style type="text/css" media="screen">
    html { overflow: auto; }


</style>
<script type="text/javascript" charset="utf-8">
    var token = '\x7b\x22link\x22\x3a\x20\x22https\x3a\x2f\x2fwww\x2edropbox\x2ecom\x2fsh\x2fq8v738nqie29asu\x2fzfZcKwpmtl\x2fstats\x5ffor\x5fprecision\x2epy\x22\x2c\x20\x22path\x22\x3a\x20\x22\x2fpython\x22\x2c\x20\x22user\x5fid\x22\x3a\x2014062712\x2c\x20\x22name\x22\x3a\x20\x22stats\x5ffor\x5fprecision\x2epy\x22\x2c\x20\x22key\x22\x3a\x20\x22q8v738nqie29asu\x22\x7d'.evalJSON();

    document.observe("dom:loaded", function () {
        scrollTo(1,0);
        scrollTo(0,0);

        DomUtil.fillVal(_('file'), 'file-or-folder');


        //Util.copy_to_clipboard_swf(token.link, 'copy_link', 'copy_done');
    });

    function copy_done () {
        Sprite.replace($("copy_link"), 'link', 'tick');
    }

    function c2d () {
        SharingModel.copy_to_dropbox('stats\x5ffor\x5fprecision\x2epy', 'q8v738nqie29asu', '\x2fstats\x5ffor\x5fprecision\x2epy', 'zfZcKwpmtl');
    }
    SharingModel.set_viewport_size();
</script>

    <script type="text/javascript">

        var _gaq = _gaq || [];
        _gaq.push(['_setAccount', 'UA-279179-2']);
        _gaq.push(['_setDomainName', document.domain]);
        _gaq.push(['_trackPageview']);

        (function() {
            var ga = document.createElement('script'); ga.type = 'text/javascript'; ga.async = true;
            ga.src = ('https:' == document.location.protocol ? 'https://ssl' : 'http://www') + '.google-analytics.com/ga.js';
            var s = document.getElementsByTagName('script')[0]; s.parentNode.insertBefore(ga, s);
        })();

    </script>
</head>

<body class="en has_preview shmodel-body" >

    

    <div id="modal-behind" style="display: none;"></div>

    <div id="modal" style="display: none;">
        <div id="modal-box">
            <a id="modal-x" onclick="javascript: Modal.hide(); Event.stop(event); return false;" href="#"></a>
            <h2 id="modal-title"></h2>
            <div id="modal-content">
            </div>
        </div>
    </div>

    <div id="modal-overlay" style="display: none;" onclick="Modal.hide(); Event.stop(event); return false;"></div>
    <div id="floaters"></div><div id="trash-can" style="display:none"></div><div id="grave-yard" style="display:none"></div>


    <div class="external-drop-indicator top" style="display: none;"></div>
    <div class="external-drop-indicator right" style="display: none;"></div>
    <div class="external-drop-indicator bottom" style="display: none;"></div>
    <div class="external-drop-indicator left" style="display: none;"></div>


    <div id="translate-div" style="display:none; margin-bottom: -10px">
        <form id="translation-suggest-form" onsubmit="TranslationSuggest.submit_suggest(event); return false;" action='/translation_suggest' method="POST">
            <input type="hidden" name="locale" value="en" />
            <input type="hidden" name="locale_url" value="https://www.dropbox.com/sh/q8v738nqie29asu/zfZcKwpmtl/stats_for_precision.py" />
            <input type="hidden" id="translation-msg-id" name="msg_id" value="" />
            <p class="clean"><label for="bad_text">Please paste or type the improper translation</label></p>
            <div id="part-one">
                <p>
                    <span class="error-message" id="bad-i18n-text-error" style="display: none">We couldn't find that string on this page.</span>
                    <textarea name="bad_text" id="bad-i18n-text" class="act_as_block textinput" rows='3' cols='40'></textarea>
                </p>
                <div id="bad-i18n-text-complete" class="autocomplete"></div>
            </div>
            <div id="part-two">
                                    <div class="emo hotbox" id="translation-msg-display" style="margin-bottom: 10px"></div>
                  <p class="clean">Original English text</p>
                  <div class="green-hotbox" id="translation-orig-msg-display" style="margin-bottom: 10px"></div>
                  <p class="clean"><label for="suggested_text">Suggested translation</label></p>
                  <p><textarea name="suggested_text" class="textinput act_as_block" rows='3' cols='40'></textarea></p>
                  <p class="clean"><label for="explanation">Explanation of the problem</label></p>
                  <p><textarea name="explanation" class="textinput act_as_block" rows='3' cols='40'></textarea></p>
                  <div class="modal-buttons">
                      <input type="button" class="freshbutton" onclick="TranslationSuggest.start_wizard(event);return false;" value="Back" id="translation-back-button"/>
                      <input type="submit" class="freshbutton-blue" onclick="TranslationSuggest.submit_suggest(event); return false;" value="Suggest translation" />
                  </div>
            </div>
        </form>
    </div>

<div id="notify-wrapper">
  <span id="notify" style="display: none;">
    <span id="notify-msg">Sorry, there was a problem loading this page.</span>
  </span>
</div>


<div id="outer-frame">
    <div id="page-header">
    </div>



    <div id="modal-locale-selector" style="display:none;">
      <ul>
        <li><a class="locale-option" data-locale="en">English</a></li>
        <li><a class="locale-option" data-locale="es">Español</a></li>
        <li><a class="locale-option" data-locale="fr">Français</a></li>
        <li><a class="locale-option" data-locale="de">Deutsch</a></li>
        <li><a class="locale-option" data-locale="ja">日本語</a></li>
      </ul>
    </div>
    <noscript><p class="center">The Dropbox website requires JavaScript.</p></noscript>
    <div id="page-content">
        <script type="text/javascript">
  document.observe('dom:loaded', function () {
    if ($('login-create-an-account')) {
      $('login-create-an-account').writeAttribute('href', '/register?signup_tag=shmodel');
    }
  });
</script>

<script type="text/javascript" charset="utf-8" src="/static/javascript/logger.js">
</script>

<script type="text/javascript">
    var shmodel_logger = new ShmodelLogger(window.location.href);
</script>

<div id="file-preview-modal" style="display:none;">
    <div class="header">
        <a class='close' href="#"><img src="&#47;static&#47;images&#47;icons&#47;icon_spacer.gif" alt="" class="sprite s_white_x " /></a>
    </div>
    <div class="preview">
        <table><tr><td class='preview-container'>
            <div class='preview-content'></div>
        </td></tr></table>
    </div>
    <div id="file-preview-menu" class="menu">
        <div class="filename">&nbsp;</div>
        <div class="actions">&nbsp;</div>
        <div class="paging">
            <a href="#" class='prev'><img src="&#47;static&#47;images&#47;icons&#47;icon_spacer.gif" alt="" class="sprite s_arrow_left_white " /></a> <span class='current_index'>&nbsp;</span> of <span class='total'>&nbsp;</span> <a href="#" class='next'><img src="&#47;static&#47;images&#47;icons&#47;icon_spacer.gif" alt="" class="sprite s_arrow_right_white " /></a>
        </div>
    </div>
    <div class="delete-file-prompt chat-bubble-bottom" style="display:none;">
        <input type="button" class="freshbutton-blue" value="Delete" id="lightbox-delete-photo"/><br/>
        <input type="button" class="freshbutton" value="Cancel" id="lightbox-delete-cancel"/>
        <div class="chat-bubble-arrow-border"></div>
        <div class="chat-bubble-arrow"></div>
    </div>
</div>


<div class="nav-header">
    <a class="logo" target="_top" href="//www.dropbox.com?src=shmodel">
        <img src="/static/images/logo-shmodel-big.png" alt="Dropbox" class="big"/>
        <img src="/static/images/logo-shmodel.png" alt="Dropbox" class="small"/>
    </a>
    <div class="filename shmodel-filename">
        stats_for_precision.py
    </div>
        <div class="buttons">
            <a id="download_button_link" href="https://dl.dropbox.com/sh/q8v738nqie29asu/zfZcKwpmtl/stats_for_precision.py?dl=1" class="freshbutton-lightblue download-button">
              <img src="&#47;static&#47;images&#47;icons&#47;icon_spacer.gif" alt="" class="sprite s_download_arrow " />Download
            </a>
              <a href="#" class='freshbutton-lightblue' id="owner-menu-button" onclick="SharingModel.toggle_owner_menu(event);"><img src="&#47;static&#47;images&#47;icons&#47;icon_spacer.gif" alt="" class="sprite s_cog " /></a>
              <div id="owner-menu" class="chat-bubble" style="display:none;">
                  <ul>
                    <li>
                        <a href="/c/browse/python/stats_for_precision.py?src=shmodel">
                            <img src="&#47;static&#47;images&#47;icons&#47;icon_spacer.gif" alt="" class="sprite s_dropbox " />Show in my Dropbox
                        </a>
                    <li>
                        <a
                               href="https://www.dropbox.com/sh/q8v738nqie29asu/x9b6fqLDX4?rem=1"
                           class="remove-link">
                          <img src="&#47;static&#47;images&#47;icons&#47;icon_spacer.gif" alt="" class="sprite s_lock " />Remove link
                        </a>
                    </li>
                  </ul>
              </div>


<div id='account-header'>
  <ul class="nav">

                <li class="ui-button">
                    <a class="header-nav-link" href="#" onclick="return false;">
                      
<span id="emsnippet-1a8b6df764f4688d"></span>
<script type="text/javascript">
document.observe('dom:loaded', function () {
  $('emsnippet-1a8b6df764f4688d').innerHTML = 'john\x20hancock'.em_snippet(20, 0.750).escapeHTML();
});
</script>
 <img src="&#47;static&#47;images&#47;icons&#47;icon_spacer.gif" alt="" class="sprite s_login_arrow link-img" />
                    </a>
                    <div class="sub-nav chat-bubble">
                        <ul>
                            <li>
                                <div class="name">john hancock</div>
                                <div class="email force-no-break">
<span id="emsnippet-980e3f5052d69297"></span>
<script type="text/javascript">
document.observe('dom:loaded', function () {
  $('emsnippet-980e3f5052d69297').innerHTML = 'selenop\x40gmail\x2ecom'.em_snippet(22, 0.250).escapeHTML();
});
</script>
</div>
                            </li>
                            <li>
                                    
                                <div class="quota-text">0.9 GB of 5.25 GB used</div>
                                <div class="quota_graph_container">
                                   <div class="quota_graph_bar " style="width: 16%;"></div>
                                </div>
                            </li>
                            <li><a href="/account" target=""><img src="&#47;static&#47;images&#47;icons&#47;icon_spacer.gif" alt="" class="sprite s_settings_16 " />Settings</a></li>
                            <li><a href="/install" target=""><img src="&#47;static&#47;images&#47;icons&#47;icon_spacer.gif" alt="" class="sprite s_install " />Install</a></li>
                            <li><a href="/plans" target=""><img src="&#47;static&#47;images&#47;icons&#47;icon_spacer.gif" alt="" class="sprite s_star_16 " />Upgrade</a></li>
                            <li><a href="/logout"><img src="&#47;static&#47;images&#47;icons&#47;icon_spacer.gif" alt="" class="sprite s_signout " />Sign out</a></li>
                        </ul>
                        <div class="chat-bubble-arrow-border"></div>
                        <div class="chat-bubble-arrow"></div>
                    </div>
                </li>
            </ul>
        </div>
        </div>
</div>

<div id="shmodel-content-area">
    

    <div id="default-content">
        <img class="bigicon" src="/static/images/icons128/page_white_code.png" />
        <div class="filename shmodel-filename">stats_for_precision.py</div>
        <div class="meta">26 hrs ago &middot; 2.52 KB</div>
        <a href="https://dl.dropbox.com/sh/q8v738nqie29asu/zfZcKwpmtl/stats_for_precision.py?dl=1" class="freshbutton-blue">
          <img src="&#47;static&#47;images&#47;icons&#47;icon_spacer.gif" alt="" class="sprite s_download_arrow " />Download
        </a>
    </div>
    <div class="filename-below shmodel-filename">stats_for_precision.py</div>
    <div class="preview-box">    <script src="/static/javascript/external/syntax.js"></script>
    <link href="/static/css/syntax.css" rel="stylesheet" type="text/css" media="all" />
    <div id="code-wrapper" class="content-shadow">
        <div id='code-loading' class="center">
            <img src="/static/images/icons/ajax-loading-small.gif" alt="" class="text-img"/>Loading...
        </div>
        <pre id="code" class="brush: py;"></pre>
    </div>
    <script type="text/javascript" charset="utf-8">
        function got_text (txt) {
            if (!txt || !txt.length) {
                $('code-wrapper').remove();
                return;
            }

            $('code-loading').remove();
            txt = Util.decode_b64(txt);
            // the new markup must contain the pre-tag to work in IE.  we CANNOT set the content of the pre
            // tag directly
            // http://stackoverflow.com/questions/451486/pre-tag-loses-line-breaks-when-setting-innerhtml-in-ie
            var class_names = $('code').className;
            $("code-wrapper").innerHTML = '<pre id="code">' + txt.escapeHTML() + '</pre>';
            $('code').addClassName(class_names);

                            SyntaxHighlighter.defaults['toolbar'] = false;
                SyntaxHighlighter.defaults['quick-code'] = false;
                SyntaxHighlighter.highlight($("code"));
            window.filled_text = 1;
        }
        document.observe("dom:loaded", function () {
            Util.add_script("https://dl.dropbox.com/sh/q8v738nqie29asu/zfZcKwpmtl/stats_for_precision.py?js_callback=got_text&b64=1");
            setTimeout(function() {
                if (!window.filled_text) {
                    $(document.body).removeClassName('has_preview');
                }
            }, 10000);
        });
    </script>
</div>
    <a class="content-flag title_bubble"
       title="Flag for copyright"
       rel="nofollow"
       href="/copyright_complaint_form?ssu=https%3A//www.dropbox.com/sh/q8v738nqie29asu/zfZcKwpmtl/stats_for_precision.py">
      Flag for copyright
    </a>

</div>

<div id="twitter-login" style="display:none;">
    <form action="#" onsubmit="Twitter.login(); return false">
        <div class="center">
            <img src="/static/images/twitter_logo.jpg" alt="twitter"/>
            <p id="twitter-status" style="margin-top:1.5em;">Before you can tweet this link, we need to connect to your Twitter account.</p>
            <p style="margin: 1.5em 0 3.5em 0;">
                <a href="#" onclick="Twitter.do_auth(); return false;">
                    <img src="/static/images/twitter_signin.png" alt=""/>
                </a>
            </p>
            <div class="modal-buttons">
                I've finished signing in to Twitter
                <input type="button" class="freshbutton-blue" onclick="Twitter.onLoginSuccessCallback();" value="Tweet my message"/>
            </div>
        </div>
    </form>
</div>

<div id="facebook-auth" style="display:none;">
    <div class="center">
        <img src="/static/images/facebook_logo.jpg" alt="facebook"/>
        <p id="facebook-status" style="margin-top:1.5em;">Before you can post this link, we need to connect to your Facebook account.</p>
        <p style="margin: 1.5em 0 3.5em 0;">
            <a href="#" onclick="FacebookOAuth.do_auth(); return false;">
                <img src="/static/images/facebook_signin.png" alt=""/>
            </a>
        </p>
        <div class="modal-buttons">
            I've finished signing in to Facebook <input  style="margin-left: 0.5em;margin-bottom: -1px;" type="button" class="freshbutton-blue" onclick="Modal.vars.action();" value="Post my message"/>
        </div>
    </div>
</div>

<div id="twitter-posting" style="display:none;">
    <div class="center">
        <img src="/static/images/twitter_logo.jpg" alt="Twitter"/>
        <h3 style="margin-top: 3em;"><img src="/static/images/icons/ajax-loading-small.gif" class="h2-img" alt=""/>Tweeting</h3>
    </div>
</div>
<div id="facebook-posting" style="display:none;">
    <div class="center">
        <br /><br />
        <img src="/static/images/facebook_logo.jpg" alt="Facebook"/>
        <br /><br />
        <h3><img src="/static/images/icons/ajax-loading-small.gif" class="h2-img" alt=""/>Posting to Facebook</h3>
        <br />
    </div>
</div>
<div id="twitter-prompt" style="display:none;">
    <form action="#" onsubmit="Twitter.custom_post(); return false;">
        <table style="width:100%; margin-bottom: 1em;">
            <tr>
                <td style="width:48px;"><img src="/static/images/twitter_default_profile.png" id='twitter-profile-image' alt=""/></td>
                <td><textarea id="twitter-msg" rows="4" maxlength="140" cols="20" style="height: 38px;" class="textinput act_as_block"></textarea></td>
            </tr>
        </table>
        <div class="modal-buttons">
            <input type="submit" class="freshbutton-blue" value="Tweet on Twitter" />
            <input type="button" class="freshbutton" value="Cancel" onclick="Modal.hide();"/>
        </div>
    </form>
</div>

<div id="post-options" style="display:none;">
    <div class="hotbox dark" style="padding: 15px 10px;">
        <table style="border-collapse:collapse; font-size: 14px; color:#095db2; margin: 0 auto;">
            <tr>
                <td colspan="2">
                <td style="padding-right: 50px;">
                    <input type="checkbox" name="facebook-checkbox" id="facebook-checkbox" checked="checked"/>
                    <label for="facebook-checkbox"><img src="&#47;static&#47;images&#47;icons&#47;icon_spacer.gif" alt="" class="sprite s_facebook text-img" />Post to Facebook</label>
                </td>
                <td>
                    <input type="checkbox" name="twitter-checkbox" id="twitter-checkbox" checked="checked"/>
                    <label for="twitter-checkbox"><img src="&#47;static&#47;images&#47;icons&#47;icon_spacer.gif" alt="" class="sprite s_twitter text-img" />Tweet on Twitter</label>
                </td>
            </tr>
        </table>
    </div>
    <div class="clearfix" style="margin-top: 1.75em">
        <span id="message-status" class="error-message"></span>
        <span id="twitter-chars-left">140</span>
        <textarea id="post-message" name="post-message" rows="4" cols="20" class="act_as_block textinput"></textarea>
    </div>
    <div class="modal-buttons">
        <input type="button" class="freshbutton-blue" id="share-this-message" value="Share this message" />
    </div>
</div>

<div id="disable-token-modal" style="display:none;">
    <p>Are you sure you want to remove this link? Once the link is removed, nobody else will be able to view this <span class="file-or-folder">&nbsp;</span>.</p>
    <div class="modal-buttons">
        <input type="button" value="Remove link" onclick="SharingModel.do_remove(Modal.vars);" class="freshbutton-blue" />
        <input type="button" value="Cancel" onclick="Modal.hide();" class="freshbutton" />
    </div>
</div>

<div id="c2d-modal" style="display: none;">
    <form action="/sm/c2d" id="c2d-form" method="post">
<input type="hidden" name="t" value="daf38c04e4" />
        <p>
            This <span class="file-or-folder">&nbsp;</span> will be instantly saved to your Dropbox and downloaded to all the computers linked to your account.
        </p>
        <div class="modal-buttons">
            <input type="submit" class="freshbutton-blue" value="Add to my Dropbox" />
            <input type="button" class="freshbutton" value="Cancel" onclick="Modal.hide();" />
        </div>
    </form>
</div>

<div id="sharing-progress" style="display: none;">
    <div class="orange-hotbox">
        <img src="/static/images/icons/ajax-loader-trans.gif" alt=""/>
        Sharing...
    </div>
</div>
<div id="sharing-posted" style="display: none;">
    <div class="green-hotbox">
        <table>
            <tr>
                <td>Message shared!</td>
                <td style="width: 5px;"></td>
                <td><a target="_blank" id="view-post"></a></td>
            </tr>
        </table>
    </div>
</div>
<div id="sharing-posted-both" style="display: none;">
    <div class="green-hotbox">
        <table>
            <tbody>
                <tr>
                    <td>Link shared!</td>
                    <td style="width: 5px;"></td>
                    <td><a target="_blank" id="view-post" href="http://www.twitter.com/"><img src="&#47;static&#47;images&#47;icons&#47;icon_spacer.gif" alt="" class="sprite s_twitter " />View tweet</a></td>
                    <td style="width: 5px;"></td>
                    <td><a target="_blank" id="view-post" href="http://www.facebook.com/profile.php?v=wall"><img src="&#47;static&#47;images&#47;icons&#47;icon_spacer.gif" alt="" class="sprite s_facebook " />View post</a></td>
                </tr>
            </tbody>
        </table>
    </div>
</div>
<div id="inline-facebook-auth" style="display:none;">
    <div class="orange-hotbox">
        You need to <a id="connect-with-facebook">connect with Facebook</a> first.
    </div>
</div>
<div id="inline-twitter-auth" style="display:none;">
    <div class="orange-hotbox">
        You need to <a id="connect-with-twitter">connect with Twitter</a> first.
    </div>
</div>

<script>
    Util.live_jag('inline-facebook-auth', 'connect-with-facebook', 'onclick', "FacebookOAuth.do_auth(); return false;");
    Util.live_jag('inline-twitter-auth', 'connect-with-twitter', 'onclick', "Twitter.do_auth(); return false;");
</script>

<script type="text/javascript">
    // Attach click listeners onto links/divs we care about
    if ($('login-hover-link')) {
        $('login-hover-link').observe('click', function() { shmodel_logger.log('click_sign_in'); });
    }

    if ($('download_button_link')) {
        $('download_button_link').observe('click', function() { shmodel_logger.log('click_download'); });
    }

    if ($('add_to_my_dropbox_link')) {
        $('add_to_my_dropbox_link').observe('click', function() { shmodel_logger.log('click_add_to_my_dropbox'); });
    }

    if ($$('.shmodel-filename').length) {
        $$('.shmodel-filename').invoke('observe', 'click', shmodel_logger.log.curry('click_filename'));
    }
</script>

    </div>

</div>
<script type="text/javascript">
    document.observe("dom:loaded", function() {
        TranslationSuggest.attach_autocomplete();
    });
</script>

<!--[if IE]>
<iframe id="hashkeeper" name='hashkeeper' width="1" height="1" style="display: none" src="/blank" onload="HashKeeper.reloading=false;"></iframe>
<![endif]-->

<div id="ieconsole" style="position: absolute; top: 0; left: 0; font-family: Courier"></div>
<div id="FB_HiddenContainer"  style="position:absolute; top:-10000px;width:0px; height:0px; left: 0;"></div>
<script>
TranslationSuggest.update_i18n_messages({});
</script>



<div id="notice-container" class="clearfix" style='display:none;'>
</div>
</body>
</html>
